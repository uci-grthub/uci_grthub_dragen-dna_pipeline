# uci_grthub_dragen-dna_pipeline

DRAGEN whole-genome sequencing (WGS) analysis in germline or somatic
tumor-only mode, driven by Snakemake (`workflow/Snakefile`, config in
`config/config.yaml`) and run through pixi.
The current config and `metadata/metadata.csv` are set up for the
261118_FrumD_WGS project.

Per sample, one DRAGEN 4.4.7 run on the hg38 pangenome (hashtable v11) does:

- adapter + quality trimming, map/align, duplicate removal
- genome-wide SNV/indel calling (+ gVCF)
- CNV calling with self-normalization
- SV calling
- repeat-expansion genotyping (ExpansionHunter, DRAGEN default catalog)
- targeted callers + star-allele / PGx calling

The small-variant VCF is then annotated with Nirvana (GRCh38), QC is aggregated
with MultiQC, and a PDF pipeline report is rendered.

## Quick start

```bash
pixi run samples            # build data/samples.csv from data/FASTQ + metadata/metadata.csv
pixi run dryrun-germline    # show planned jobs (or dryrun-tumor-only)
pixi run run-germline       # run; DRAGEN jobs serialized (resources: dragen=1)
pixi run run-tumor-only     # somatic tumor-only mode
```

`mode` has no default; `pixi run dryrun` / `run` only work once `mode:` is set
in `config/config.yaml`.

## Samples

FASTQs live in `data/FASTQ` (`<run>-<lane>-...-<i7>-<i5>-R{1,2}.fastq.gz`).
`src/make-sample-file.sh --metadata metadata/metadata.csv` names each pair by
matching the `<i7>-<i5 NovaSeq v1.5>` barcodes in the filename against
`metadata/metadata.csv`:

| sample | i7 | i5 (NovaSeq v1.5) |
|---|---|---|
| PT2 | ACCACACG | CAACCACG |
| PT5 | ACCAGGAC | CCGCTACA |

Without `--metadata` the sample name is the filename prefix before `-R1`.
Unmatched FASTQs are skipped with a warning.

## Config toggles (`config/config.yaml`)

| key | default | DRAGEN flags |
|---|---|---|
| `enable_cnv` | `true` | `--enable-cnv true --cnv-enable-self-normalization true` |
| `enable_sv` | `true` | `--enable-sv true` |
| `enable_repeats` | `true` | `--repeat-genotype-enable true` |
| `enable_pgx` | `true` | `--enable-targeted true --enable-star-allele true --enable-pgx true` |
| `emit_gvcf` | `true` | `--vc-emit-ref-confidence GVCF --vc-enable-vcf-output true` |
| `vc_filter` | `""` | `--vc-hard-filter` (omitted when empty) |
| `enable_panel` | `false` | gene-panel extraction (not a DRAGEN flag; see below) |

## Gene panel (optional, `enable_panel`)

Panel-agnostic gene-panel extraction (machinery inherited from the KainZ exome
project). When `enable_panel: true`,
genes in `panel_csv` (col 1 `gene_symbol`) are mapped to coordinates via the
GENCODE GTF. Each annotated VCF is then subset to those genes, and the subsets
are merged across samples. The merged sites are annotated once with Nirvana
JSON (gnomAD, ClinVar, dbSNP, REVEL, phyloP) and loaded into:

- `results/<panel_name>_variants.sqlite` (+ `.xlsx`, `.vcf.gz`)
- `results/<panel_name>_report.pdf`

The panel is defined by three config keys:

| key | meaning |
|---|---|
| `panel_csv` | gene list CSV (header row; `gene_symbol` = column 1, the only required column) |
| `panel_name` | basename for `results/<panel_name>_*`, shown in reports |
| `panel_description` | paragraph printed in the PDF report + xlsx `Overview` sheet |

Optional CSV columns (any order after column 1, missing = empty): `ensembl_gene_id`,
`category`, `priority`, `mechanism`, `key_variants` (rsIDs here set
`is_key_variant`), `evidence_note`, `exome_capture_note`.

To use a different panel, point `panel_csv` at the new CSV and change
`panel_name` and `panel_description`. Each `panel_name` writes its own
`results/<panel_name>_*` files, so the old panel's results are kept.

Current default: `data/pain_gene_panel_from_zeev_kain.csv` (pain panel copied
from the KainZ project). `data/` is gitignored, so copy the file in again on a
fresh clone.
`pixi run pain-genes` rebuilds the GO + CPIC/PharmGKB sourced list in
`data/pain_genes.tsv`, which has a different format and is not read by the workflow.

One-off run without editing config: `pixi run snakemake --cores all --resources dragen=1 --config enable_panel=True`.

## Outputs

`mode` (required: `germline` or `tumor_only`) picks the output root. Each
root holds its own `logs/`, `benchmarks/`, `multiqc/` and `versions.yaml`.

```
output/germline/<sample>/            mode: germline
output/somatic/align/<sample>/       mode: tumor_only, shared tumor BAM
output/somatic/<calib>/<sample>/     mode: tumor_only, calls (calib: default, auto)
results/germline/, results/somatic/  reports (+ panel)
```

Germline, per sample under `output/germline/<sample>/`:

- `<sample>.bam` (+ `.bai`) — mapped, deduplicated reads
- `<sample>.hard-filtered.vcf.gz` — small variants; `.hard-filtered.gvcf.gz` gVCF
- `<sample>.annotated.vcf.gz` — Nirvana-annotated small variants
- `<sample>.cnv.vcf.gz`, `<sample>.sv.vcf.gz`, `<sample>.repeats.vcf.gz`
- star-allele / PGx / targeted-caller outputs, DRAGEN QC metrics

Tumor-only, per sample under `output/somatic/<calib>/<sample>/`:

- `<sample>.somatic_params.yaml` — callability threshold / CNV bin width used
- `<sample>.hard-filtered.vcf.gz` — somatic small variants (INFO `GermlineStatus`)
- `<sample>.cnv.vcf.gz`, `<sample>.sv.vcf.gz`, TMB / HRD metrics,
  callable-region BEDs

See `docs/germline-vs-somatic.md` for how the two modes and the two
calibration variants differ.

Cohort-level (per mode root):

- `<root>/multiqc/multiqc_report.html` — aggregated QC
- `<root>/versions.yaml` — live tool/reference versions
- `results/germline/pipeline_report.pdf` — methods + versions + outputs

## Development

Format the Snakefile with `pixi run snakefmt workflow/Snakefile`.

Legacy standalone drivers (`src/run-nirvana.sh`,
`src/WGS_bwa_gatk.sh`, `src/run_chen.sh`) are not used by the workflow.
See `docs/updating-nirvana-cache.md` for refreshing Nirvana annotation data.
