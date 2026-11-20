# Deliverables: PI3K + RAS/MAPK mutation profile, endometrial cancer organoids (PT2, PT5)

## Context
Client (FrumD lab) has a manuscript under review. Reviewers want the mutation status of
oncogenes and tumor suppressors in the PI3K and RAS/MAPK pathways for their patient-derived
endometrial cancer organoids. Patients are de-identified, no clinical sequencing, **no matched
normal**. We have WGS of 2 organoids (PT2, PT5; ~31–34x; both female, Y cov 0).

Current pipeline state: DRAGEN 4.4.7 **germline** mode, Nirvana annotation, MultiQC,
pipeline report. Germline mode is wrong tool for tumor: assumes diploid 0/0.5/1 VAF, misses
subclonal and LOH-skewed somatic calls. Need DRAGEN **tumor-only somatic** mode, then a
pathway gene panel using existing panel machinery (`enable_panel`, `src/panel_tables.py`,
`src/export_panel_xlsx.py`, `src/make_panel_report.sh`).

## Proposed deliverables (to client)

1. **Supplementary Table (xlsx), reviewer-ready** — one row per variant in panel genes, per sample:
   gene, pathway (PI3K / RAS-MAPK), role (oncogene / TSG), HGVS c./p., consequence, VAF, depth,
   gnomAD AF, ClinVar, COSMIC/OncoKB hotspot flag, tier (AMP/ASCO/CAP: Tier I/II likely
   driver, Tier III VUS), germline-vs-somatic likelihood note. Sheets: `Overview`/README,
   `Drivers` (filtered), `All_panel_variants`, `Gene_coverage`, `Methods`.
2. **Oncoprint figure** (PDF + PNG, publication-ready) — genes x {PT2, PT5}, colored by
   alteration type (missense hotspot, truncating, amplification, deep deletion/LOH),
   grouped by pathway. Possibly a pathway schematic with hits overlaid.
3. **Copy-number summary for panel genes** — PTEN loss, PIK3CA/KRAS/ERBB2 amplification,
   NF1 loss etc. from DRAGEN somatic CNV (tumor-only) — fold into table + oncoprint.
4. **Molecular-subtype context (ProMisE/TCGA-style)** — POLE exonuclease domain hotspots,
   MSI status (DRAGEN `--msi-command tumor-only`), TMB, TP53 status. Reviewers of endometrial
   papers commonly expect this; cheap to add from same run. Flag as optional.
5. **Methods paragraph + limitations** (manuscript-ready text) — DRAGEN version, reference,
   tumor-only calling, filtering (gnomAD AF < 0.1%, PON), annotation sources, coverage, and
   explicit caveat that germline vs somatic cannot be separated without matched normal.
6. **Per-gene coverage QC** — mean depth + % bases ≥20x over coding exons of each panel gene,
   so absent calls can be stated as "no mutation detected at ≥Nx".
7. **Data package** — somatic VCFs (annotated), BAMs (on request), MultiQC report,
   `versions.yaml`, plus data-availability guidance (dbGaP/SRA controlled access; ask
   whether consent allows deposit).

## Default gene panel (`data/pi3k_rasmapk_panel.csv`)
- PI3K: PIK3CA, PIK3R1, PIK3R2, PTEN, AKT1, AKT2, AKT3, MTOR, TSC1, TSC2, PIK3CB, INPP4B, STK11, RICTOR, RPTOR
- RAS/MAPK: KRAS, NRAS, HRAS, BRAF, RAF1, ARAF, MAP2K1, MAP2K2, MAPK1, NF1, SOS1, PTPN11, RASA1, SPRY2, DUSP4
- Upstream RTK: ERBB2, ERBB3, EGFR, FGFR2, FGFR3
- Subtype context: POLE, TP53, ARID1A, CTNNB1, MLH1/MSH2/MSH6/PMS2
(columns: gene_symbol, category=pathway, mechanism=oncogene|TSG, key_variants=hotspots)

## Implementation (repo changes)
- `config/config.yaml`: add `mode: somatic_tumor_only` toggle; panel_csv/name/description
  for new panel; `enable_panel: true`.
- `workflow/Snakefile`: DRAGEN tumor-only rule (`dragen_somatic`, `mode: tumor_only`)
  (`--tumor-fastq-list` / `--tumor-bam-input` re-use existing BAMs to skip realignment,
  `--enable-variant-caller true`, germline tagging + population filter via
  `--vc-enable-germline-tagging`/`--vc-systematic-noise` PON, `--enable-cnv` somatic
  tumor-only, `--msi-command tumor-only`, `--enable-tmb true`). Confirm flags via
  `dragen-docs` skill before coding.
- Nirvana on somatic VCF (existing `src/run-nirvana.sh` rule reused).
- Panel extraction reuses existing machinery; extend `src/panel_tables.py` with VAF,
  hotspot flag, tier columns; extend `src/export_panel_xlsx.py` with `Drivers` and
  `Gene_coverage` sheets.
- New `src/oncoprint.R` or python (matplotlib) rule for figure.
- Format with `snakefmt`; envs via pixi (no venvs).

## Verification
- `pixi run dryrun` shows new somatic + panel jobs only.
- Positive controls: check known endometrial hotspots recovered if present (PIK3CA
  E545K/H1047R, KRAS G12/G13, PTEN truncations); manually inspect each driver call in IGV
  (BAM + somatic VCF), screenshots into report appendix.
- Cross-check germline VCF: driver also present at ~50% VAF + gnomAD common → flag germline.
- Per-gene coverage ≥20x over coding exons; report any gaps.

## Open questions for client
- Gene list OK, or reviewers named specific genes?
- Want subtype context (POLE/MSI/TMB/TP53) included?
- Any parental tumor tissue or normal (blood/adjacent) available? Would resolve germline/somatic.
- Organoid passage number (culture drift caveat).
- Data deposit plans / consent scope.
