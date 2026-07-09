# uci_grthub_dragen

DRAGEN drivers for running multiple samples through mapping, variant calling,
and Nirvana annotation. Two pipelines:

- **Genome (WGS)** — `src/dragen-run.sh`
- **Exome** — `src/dragen-run-exome.sh` (target-restricted, padded, on-target QC)

Both read a `sample,R1,R2` CSV. `src/make-sample-file.sh` builds that CSV by
pairing R1/R2 FASTQs in a directory (sample name = filename up to the `-R1`
tag; unpaired files skipped with a warning; header `#`-commented so the runners
skip it):

```bash
./src/make-sample-file.sh --abs --output data/samples.csv
```

Options: `--input-dir` (default `data`), `--output` (default stdout), `--abs`,
`--r1-tag` / `--r2-tag` (defaults `R1` / `R2`).

---

## Genome (WGS) — `dragen-run.sh`

Whole-genome mapping + variant calling over the full reference. No target BED.

- **Reference**: hg19 hash table
  (`/staging/human/reference/hg19/hg19.fa.k_21.f_16.m_149`)
- Adapter + quality trimming, mark/remove duplicates
- Nirvana annotation (GRCh37) of the hard-filtered VCF

```bash
./src/dragen-run.sh \
  --reference /staging/human/reference/hg19/hg19.fa.k_21.f_16.m_149 \
  --input-dir /staging/Dragen_Dev/data/FASTQ \
  --sample-file samples.csv \
  --output-base /staging/Dragen_Dev/results
```

Dry-run first:

```bash
./src/dragen-run.sh --reference /path/ref --input-dir /data/fastq \
  --sample-file samples.csv --output-base /tmp/out --dry-run
```

Key options: `--reference`, `--input-dir`, `--sample-file`, `--output-base`,
`--rgid`, `--vc-filter`, `--enable-variant-caller true|false`, `--dry-run`,
`--create-only`.

---

## Exome — `dragen-run-exome.sh`

Restricts variant calling to the IDT xGen Exome Research Panel v2 target BED
(hg38), pads targets, reports on-target coverage QC, and writes BAMs.

- **Reference**: hg38 hash table (`/staging/human/reference/hg38/hg38_hash_table`)
- **BEDs**: `/staging/human/reference/hg38/exome_beds/idt_x_gen_exome_research_panel_v2`
  - kit: `xgen-exome-hyb-panel-v2-{targets,probes}-hg38.bed`
  - spike-in: `spike-in-xgen-exome-v2-{targets,probes}-hg38.bed`
- Adapter + quality trimming, mark/remove duplicates
- Variant calling over target BED, padded 100 bp (default)
- On-target coverage QC (target BED = region 1, probe BED = region 2)
- Mapped reads written to `<sample>.bam` (+ `.bai`)
- Nirvana annotation (GRCh38) of the hard-filtered VCF

Kit only:

```bash
./src/dragen-run-exome.sh \
  --sample-file data/samples.csv \
  --output-base results
```

With spike-in panel (this prep used it — merges kit + spike-in BEDs into a
union for both calling and QC; requires `bedtools`, installed in the miniforge
base env). Merged BEDs written to `<output-base>/combined-{targets,probes}-hg38.bed`:

```bash
./src/dragen-run-exome.sh \
  --sample-file data/samples.csv \
  --output-base results \
  --with-spike-in
```

### Exome options

| Option | Default | Meaning |
|--------|---------|---------|
| `--reference PATH` | hg38 hash table | DRAGEN hash-table reference dir |
| `--input-dir PATH` | (none) | FASTQ dir; used when sample file has basenames |
| `--sample-file PATH` | **required** | CSV `sample,R1,R2` (`#` lines ignored) |
| `--output-base PATH` | **required** | Base output dir (per-sample subdir created) |
| `--target-bed PATH` | xGen v2 targets | Target BED for variant calling |
| `--target-padding N` | `100` | Bases padded each side of targets (calling only) |
| `--probe-bed PATH` | xGen v2 probes | Probe BED for coverage QC (region 2) |
| `--with-spike-in` | off | Merge spike-in BEDs into kit BEDs (union) |
| `--spike-target-bed PATH` | spike-in targets | Spike-in target BED |
| `--spike-probe-bed PATH` | spike-in probes | Spike-in probe BED |
| `--rgid STR` | `Illumina_RGID` | Read group ID |
| `--vc-filter STR` | *empty* | Hard-filter string; omitted when empty |
| `--enable-variant-caller true\|false` | `true` | Toggle variant calling |
| `--dry-run` | off | Print commands, do not execute |
| `--create-only` | off | Create dirs + print commands, do not run DRAGEN |

### Exome outputs (per sample, under `<output-base>/<sample>/`)

- `<sample>.bam` + `<sample>.bam.bai` — mapped, dedup'd reads
- `<sample>.hard-filtered.vcf.gz` — variants
- coverage QC reports (`*.qc-coverage-region-*`)
- Nirvana annotation files (if run)

### Notes

- Padding applies to variant-calling regions only, not the QC BED.
- `--vc-filter` empty by default → no `--vc-hard-filter` passed. Provide a value
  to enable DRAGEN hard filtering, e.g.
  `--vc-filter 'SNP filter:snp:QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0; INDEL filter:indel:QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'`
- `--with-spike-in` needs `bedtools` on `PATH`.
- Nirvana paths default to GRCh38; override via `NIRVANA_*` environment variables.
