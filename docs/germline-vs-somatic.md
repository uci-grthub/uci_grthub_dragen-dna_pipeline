# Germline vs tumor-only somatic DRAGEN invocation

Comparison of the two DRAGEN 4.4.7 rules in `workflow/Snakefile`:
`dragen` (`mode: germline`) and `dragen_somatic` (`mode: tumor_only`, DRAGEN
"DNA Somatic Tumor-Only Solid WGS" recipe).

**Shared:** adapter + quality trimming (Q15, min length 20), FASTQ input,
read group, map/align, BAM output, duplicate handling.

## Reference

| | Germline | Somatic |
|---|---|---|
| Hashtable | `pangenome_reference` (graph hg38, v11, CNV-enabled) | linear hg38. Was `hg38_hash_table_v11` (no `kmer_cnv.bin`, so the first run failed); moving to Illumina's `hg38-alt_masked.cnv.hla.methyl_cg.methylated_combined.rna-11-r5.0-1` |

DRAGEN recommends a linear hashtable for somatic calling. The somatic BAM is
therefore not the same alignment as the germline BAM.

## Input and duplicates

| | Germline | Somatic |
|---|---|---|
| Reads | `-1/-2`, `--RGSM` | `--tumor-fastq1/2`, `--RGSM-tumor` (tumor input with no normal is what selects tumor-only mode) |
| Duplicates | `--remove-duplicates true` (removed from BAM) | `--enable-duplicate-marking true` (flagged, kept in BAM) |

Callers ignore marked duplicates, so calls are not affected by this difference.

## Small-variant calling

This is the part that matters for the PI3K / RAS-MAPK results.

**Germline**
- Diploid genotyping: each variant must fit 0/1 or 1/1, so variants are
  expected near 50% or 100% of reads.
- gVCF output; no hard filter.

**Somatic**
- No ploidy model. Calls down to 3% VAF (`--vc-target-vaf 0.03`), so
  subclonal variants and variants with allele fractions skewed by copy loss
  are kept.
- `--vc-systematic-noise`: removes recurrent sequencing / alignment
  artifacts using Illumina's WGS hg38 v2.0.0 baseline (40 fresh-frozen
  samples, 70x).
- `--vc-enable-germline-tagging` + `--variant-annotation-data`: tags
  variants common in gnomAD / 1000 Genomes as likely germline, standing in
  for the missing matched normal.
- Output VCF carries FILTER status and per-variant allele fractions; no
  germline-style genotypes, no gVCF.
- Germline tagging marks variants with INFO `GermlineStatus`
  (`Germline_DB`, `Germline_Proxi`, `Somatic`, `Somatic_Hotspot`, ...); it
  does not change FILTER. Somatic candidates = PASS **and**
  `GermlineStatus=Somatic*`.
- Germline-tagging caveat for organoids: they are ~100% tumor, so clonal
  heterozygous somatic variants sit near 50% VAF like germline ones. The TMB
  proxy filter (`--tmb-enable-proxi-filter`) tags such calls
  `Germline_Proxi` (24k PASS calls in PT5). The panel step must report
  `Germline_Proxi` calls in pathway genes as "somatic or germline,
  unresolved", not drop them.

### Two somatic calibration variants

The somatic run is two DRAGEN passes: `dragen_somatic_align` (map/align,
coverage metrics; shared, `output/somatic/align/<sample>/`) then
`dragen_somatic` (calling from that BAM), once per entry in
`somatic_calibrations`, into `output/somatic/<calib>/<sample>/`. In between,
`somatic_params` writes `<sample>.somatic_params.yaml` with the flags used:

| | `default` | `auto` |
|---|---|---|
| `--vc-callability-tumor-thresh` | DRAGEN default (50 reads) | clamp(round(0.5 × median autosomal cov), 10, 50); PT5 13, PT2 14 |
| `--cnv-interval-width` | DRAGEN default (1000 bp) | clamp(ceil_1kb(1000 × 100 / median), 1000, 10000); PT5 and PT2 4000 bp |
| `--vc-germline-tag-hotspots` | DRAGEN default | `false` (COSMIC hotspots never tagged germline) |

Why: DRAGEN's defaults assume ~100x tumors.
- Callability 50 on PT5 (median 25x) left 5 Mb of genome and 0.04 Mb of the
  TMB region callable, so TMB came out as 0. The threshold does **not**
  change calls: PT2 and PT5 VCFs at 50 vs 13/14 have identical variant
  counts. It changes only callable-region BEDs and TMB.
- 1 kb CNV bins at 25x: `WARNING: Autosomal median is too low. Consider
  increasing --cnv-interval-width`, then the model plausibility check failed
  (`ModelCallsTooManyHomdels`) and the CNV VCF had 0 records (purity NA,
  ploidy 0, no HRD).

The `default` results were produced by the first one-pass run (same reads,
same reference, DRAGEN defaults) and moved into `output/somatic/default/`.

## CNV

| | Germline | Somatic |
|---|---|---|
| Flags | `--enable-cnv --cnv-enable-self-normalization` | same, plus `--cnv-population-b-allele-vcf` (1000G phase 1 high-confidence SNPs) |
| Model | integer copy states, diploid genome assumed | allele balance at population SNP sites gives LOH, purity and ploidy without a normal; somatic segmentation (`aslm`); reports copy-neutral LOH |

The somatic CNV calls are what the PTEN-loss and amplification
(PIK3CA / KRAS / ERBB2) deliverable needs.

## SV

| | Germline | Somatic |
|---|---|---|
| Flags | `--enable-sv` (germline mode) | `--enable-sv` (tumor mode) + `--sv-systematic-noise` (WGS hg38 v3.1.0) |
| Model | germline SV genotyping | somatic scoring, low allele fractions allowed |

## Somatic-only outputs

- `--enable-hrd`: homologous-recombination-deficiency score, built from the CNV calls.
- `--enable-tmb --tmb-enable-proxi-filter`: tumor mutational burden with a
  VAF-based germline proxy filter. Tumor-only TMB still tends to run high
  because some germline variants get through.
- MSI: wired in but off (`enable_msi: false`). The only WGS normal baseline
  Illumina ships is FFPE-derived, and ~31-34x coverage leaves few sites above
  threshold.
- `--enable-metrics-json`.

## Germline-only outputs

Dropped from the somatic run because they mean nothing in tumor-only
analysis:

- repeat-expansion genotyping (ExpansionHunter)
- targeted callers, star-allele and PGx calls
- gVCF

## Practical consequence

Germline variants in the somatic VCF are **tagged, not removed**. A
tumor-only call near 50% VAF that is absent from gnomAD may still be a rare
germline variant. The germline VCF of the same organoid cannot separate
the two, because both runs use tumor DNA. Only a matched normal (blood or
adjacent tissue) can. The manuscript methods paragraph needs to state this.
