# Variant counts: germline vs somatic default vs somatic auto

PT2 and PT5 endometrial cancer organoids, WGS (~25-27x median autosomal
coverage), DRAGEN 4.4.7. Generated 2026-11-20 from:

| Run | Files | Reference | Settings |
|---|---|---|---|
| germline | `output/germline/<s>/<s>.hard-filtered.vcf.gz` | hg38 pangenome (graph) | germline WGS |
| somatic default | `output/somatic/default/<s>/` | hg38 alt-masked linear (CNV-enabled) | tumor-only, DRAGEN defaults: callability 50, CNV bins 1 kb |
| somatic auto | `output/somatic/auto/<s>/` | same | tumor-only; callability PT2 14 / PT5 13, CNV bins 4 kb, hotspots never germline-tagged |

The two somatic runs call from the same reads and alignment. See
`docs/germline-vs-somatic.md` for what each setting does.

## Small variants (genome-wide)

| Run | Sample | Total | PASS | PASS SNV | PASS indel |
|---|---|---|---|---|---|
| germline | PT2 | 5,128,109 | 5,094,098 | 4,131,048 | 968,264 |
| somatic default | PT2 | 5,003,049 | 4,452,336 | 3,640,871 | 736,793 |
| somatic auto | PT2 | 5,003,049 | 4,452,336 | 3,640,871 | 736,793 |
| germline | PT5 | 5,099,481 | 5,067,534 | 4,110,550 | 961,982 |
| somatic default | PT5 | 4,962,745 | 4,418,246 | 3,615,970 | 727,948 |
| somatic auto | PT5 | 4,962,745 | 4,418,246 | 3,615,970 | 727,948 |

- The somatic caller passes ~13% fewer variants than germline (stricter
  filtering: systematic noise, weak evidence, too few supporting reads).
- Default and auto calls are identical record for record: the callability
  threshold and CNV bin width do not change small-variant calls.
- Germline calls are not removed in tumor-only mode, only tagged (next table).

## PASS calls by germline tag (somatic runs)

| Run | Sample | Germline_DB | Germline_Proxi | Somatic | Somatic_Hotspot | untagged |
|---|---|---|---|---|---|---|
| default | PT2 | 4,388,278 | 45,085 | 5,921 | 0 | 13,052 |
| auto | PT2 | 4,387,902 | 45,083 | 5,919 | 380 | 13,052 |
| default | PT5 | 4,379,912 | 24,475 | 1,753 | 0 | 12,106 |
| auto | PT5 | 4,379,565 | 24,473 | 1,753 | 349 | 12,106 |

- `Germline_DB`: in gnomAD / 1000 Genomes. `Germline_Proxi`: judged germline
  from VAF (unreliable in ~100% tumor organoids, where clonal somatic
  variants also sit near 50% VAF). Untagged: mostly MNV / complex records.
- The only default vs auto difference is `--vc-germline-tag-hotspots false`:
  ~350-380 COSMIC hotspot positions per sample that default tagged
  `Germline_DB` (or `Somatic`) are `Somatic_Hotspot` in auto. Most are
  population SNPs that happen to sit on a hotspot position; the pathway
  filter below removes them by gnomAD AF.

## Callable region, TMB, CNV

| Run | Sample | Callable genome | TMB region | TMB (mut/Mb) | CNV gains / losses | Passing CNVs | Ploidy | Purity |
|---|---|---|---|---|---|---|---|---|
| default | PT2 | 7 Mb | 0.06 Mb | 0.00 | 0 / 0 | 0 | 0.00 | NA |
| auto | PT2 | 2,779 Mb | 34.26 Mb | 1.37 | 86 / 81 | 0 | 2.05 | NA |
| default | PT5 | 5 Mb | 0.04 Mb | 0.00 | 0 / 0 | 0 | 0.00 | NA |
| auto | PT5 | 2,765 Mb | 33.97 Mb | 0.24 | 68 / 90 | 0 | 2.00 | NA |

- Default TMB and callable regions are unusable at ~25x (50-read threshold).
- Default CNV aborted (`Autosomal median is too low`, 0 records). Auto CNV
  segments the genome, but no event passes DRAGEN's filters and purity is not
  estimated. No segment > 1 Mb deviates > 10% from a copy ratio of 1.0 in
  either sample.

## PI3K / RAS-MAPK panel (auto, `results/somatic/pi3k_rasmapk_*`)

43-gene panel (`config/panels/pi3k_rasmapk_panel.csv`). PASS calls in panel
genes, distinct variants per sample:

| Sample | likely germline (DB) | unresolved (Germline_Proxi) | somatic |
|---|---|---|---|
| PT2 | 6,289 | 98 | 31 |
| PT5 | 6,367 | 39 | 21 |

Protein-altering (high / moderate impact, canonical transcript), PASS,
gnomAD AF < 0.1%, not DB-germline (`pathway_hits` view):

| Gene | Sample | Protein change | VAF | Depth | Tag | ClinVar |
|---|---|---|---|---|---|---|
| PIK3CA | PT2 | p.Arg88Gln | 0.13 | 32 | Somatic_Hotspot | Pathogenic |
| PIK3CA | PT2 | p.Ile112_Leu113dup | 0.15 | 34 | Somatic_Hotspot | — |
| KRAS | PT5 | p.Gly12Ala | 0.39 | 18 | Somatic_Hotspot | Pathogenic; Pathogenic/Likely pathogenic |
| RASA1 | PT2 | p.Pro136Leu | 0.50 | 22 | Germline_Proxi (unresolved) | Conflicting; VUS (gnomAD AF 3.3e-05) |

Only one protein-altering, rare panel call was filtered out by DRAGEN: BRAF
p.Tyr472Cys in PT2 (VAF 0.04, depth 23, `too_few_supporting_reads`). No PTEN
call passes or fails these criteria in either sample.

Copy ratio over panel genes (CNV segment means, auto): PTEN 1.00 / 1.01,
PIK3CA 1.06 / 1.00, KRAS 1.00 / 1.00, NF1 1.00 / 1.00, ERBB2 0.98 / 0.99
(PT2 / PT5). All `segmentMean`-filtered, i.e. not distinguishable from
copy-neutral.

## Caveats

- No matched normal: somatic vs germline is inferred, not proven. RASA1
  p.Pro136Leu at VAF 0.50 may be a rare germline variant.
- PIK3CA calls in PT2 at VAF 0.13-0.15 are subclonal or the organoid is not
  ~100% tumor; verify in IGV (`output/somatic/align/PT2/PT2_tumor.bam`).
- ~25x coverage: variants below ~15% VAF have <= 4 supporting reads and are
  under-called.
