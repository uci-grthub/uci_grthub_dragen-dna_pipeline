# Updating the Nirvana / ICA annotation cache

The bundled Illumina Connected Annotations (ICA) binary is **3.26.0** (under
`/opt/dragen/4.4.7/share/nirvana/`).

**Status: resolved (2026-07-09).** Transcript annotation now works. The pipeline
points at the ICA-3.26 cache and a de-duplicated SD dir (paths below). This doc
records what was wrong, how it was fixed, and how to repeat the update when the
data is next refreshed.

## Symptom (what was broken)

The old on-disk cache (`.../Cache/GRCh38/Both.transcripts.ndb`, dated Oct 2025)
is an **old-format cache that ICA 3.26 silently rejects** — it "loads" instantly
but contributes no transcript data source, so annotated output had **no gene /
HGVS / consequence** fields. SD-dir sources (gnomAD, ClinVar, dbSNP, SpliceAI,
REVEL …) still worked; only transcripts were missing.

## How to tell the cache is broken

Run ICA on any VCF and look at the **Data Source** table it prints:

```bash
N=/opt/dragen/4.4.7/share/nirvana/Nirvana
$N -c <cache-dir> -r <ref.dat> --sd <sd-dir> -i some.vcf.gz -o /tmp/probe
```

- **Broken:** no `GeneModels` rows for `Ensembl` / `RefSeq` / `HGNC`.
  `zcat /tmp/probe.json.gz | grep -c '"transcripts"'` = 0.
- **Working:** table shows `Ensembl GeneModels Loaded`, `RefSeq GeneModels
  Loaded`, `HGNC GeneModels Loaded`, `GenomeAssembly Loaded`, and variants carry
  a `transcripts` array (annotated VCF records carry `CSQ=`).

## Current working configuration

`config/config.yaml` (and the `NIRVANA_*` defaults in `src/run-nirvana.sh` /
`src/dragen-run-exome.sh`):

```yaml
nirvana_cache: "/staging/human/reference/hg38/annotation_data_files/Cache"
nirvana_dat:   "/staging/human/reference/hg38/annotation_data_files/References/Homo_sapiens.GRCh38.Nirvana.dat"
nirvana_sd:    "/staging/human/reference/hg38/annotation_data_files/SupplementaryAnnotation/GRCh38_ica326"
```

Two non-obvious points, both learned the hard way (see below):

- **`nirvana_cache` is the parent `Cache/`, NOT `Cache/GRCh38/`.** The ICA-3.26
  transcript cache is a set of files in the `Cache/` root — `GRCh38.Ensembl.ndb`,
  `GRCh38.RefSeq.ndb`, `GeneSymbols.ndb` (+ `.idx`). The legacy
  `Cache/GRCh38/Both.transcripts.ndb` still sits in the subdir; pointing at the
  subdir gets you the rejected old cache.
- **`nirvana_sd` is `GRCh38_ica326`, a de-duplicated copy of `GRCh38/`.** See the
  duplicate-key section.

## Steps to refresh the data (repeat this when updating)

Credentials live at `~/.ilmnAnnotations/credentials.json` (tied to the Illumina
ICA account/license; keep `chmod 600`, never commit). The download license is
written to `~/.ilmnAnnotations/premium.lic`.

1. **Generate a fresh download config** — do NOT reuse an old one:

   ```bash
   DM=/opt/dragen/4.4.7/share/nirvana/DataManager
   $DM make-config -r GRCh38 -o ~/.ilmnAnnotations
   # writes ~/.ilmnAnnotations/GRCh38_annotation_config.json
   ```

   **Gotcha:** a stale/hand-trimmed config can be missing the `GeneModels`
   entries (`Ensembl`, `RefSeq`, `HGNC`) and `GenomeAssembly` — those *are* the
   transcript cache. Downloading with such a config silently reproduces the
   broken state. Confirm they are present:

   ```bash
   grep -c GeneModels ~/.ilmnAnnotations/GRCh38_annotation_config.json   # expect 3
   ```

2. **Download** (large, ~28 GB; needs network):

   ```bash
   $DM download -r GRCh38 \
     --credentials-file ~/.ilmnAnnotations/credentials.json \
     --versions-config  ~/.ilmnAnnotations/GRCh38_annotation_config.json \
     --dir /staging/human/reference/hg38/annotation_data_files \
     --thread 4
   ```

   Downloading in place upgrades the transcript cache (`Cache/GRCh38.*.ndb`) and
   `References/…Nirvana.dat`, and writes new SA versions into
   `SupplementaryAnnotation/GRCh38/`.

3. **Rebuild the de-duplicated SD dir.** `DataManager` *adds* new SA versions but
   does **not** remove the superseded old ones, so `SupplementaryAnnotation/GRCh38`
   ends up with two versions of some sources (e.g. `ClinVar_20240630` +
   `ClinVar_20260328`, `dbSNP_156` + `dbSNP_157`, and old `.nsa` twins of sources
   that moved to `.esa`: COSMIC, SpliceAi, PrimateAI-3D). ICA then aborts with:

   ```
   ERROR: Duplicate variant-level (.nsa) JSON keys found for: <source>
   ```

   Fix by pointing `--sd` at a dir that has exactly one version per source.
   `GRCh38_ica326` is a symlink farm holding the newest of each (kept
   `dbSNP_151_globalMinor` — that's the distinct RefMinor source, not a dbSNP
   dup). To regenerate it, keep the newest file per source key and drop the
   superseded ones; verify no source has two versions:

   ```bash
   SD=/staging/human/reference/hg38/annotation_data_files/SupplementaryAnnotation
   ls "$SD"/GRCh38_ica326/*.nsa "$SD"/GRCh38_ica326/*.esa \
     | xargs -n1 basename | sed -E 's/_[0-9].*//' | sort | uniq -d
   # (only 'dbSNP' should appear — globalMinor vs 157, which are different keys)
   ```

4. **Verify** transcripts load and populate:

   ```bash
   N=/opt/dragen/4.4.7/share/nirvana/Nirvana
   D=/staging/human/reference/hg38/annotation_data_files
   $N -c "$D/Cache" \
      -r "$D/References/Homo_sapiens.GRCh38.Nirvana.dat" \
      --sd "$D/SupplementaryAnnotation/GRCh38_ica326" \
      --output-format json \
      -i output/<sample>/<sample>.hard-filtered.vcf.gz -o /tmp/probe
   zcat /tmp/probe.json.gz | grep -c '"transcripts"'   # should be > 0
   ```

   Note: a working run is much slower than a broken one (~2–3 min/sample vs a
   couple of seconds) because it now does real transcript annotation.

5. **Re-annotate** the batch:

   ```bash
   pixi run run                                              # via Snakemake
   # or
   bash src/run-nirvana.sh --output-base output --force      # standalone
   ```

   DRAGEN is skipped (BAMs/VCFs already exist); only the `nirvana` rule re-runs.

## Notes

- `credentials.json` is a secret — keep it `chmod 600`, never commit it.
- The old `GRCh38_curated` SD workaround (dropped DANN / GERP / PhyloP-primate to
  avoid missing-`.idx` errors) is **no longer used** — the fresh download ships
  complete `.idx` sidecars. `nirvana_sd` now points at `GRCh38_ica326`.
- Basic-tier sources (gnomAD, ClinVar, dbSNP, GERP, DANN, phyloP, TOPMed,
  1000 Genomes, DECIPHER, FusionCatcher, MITOMAP, MultiZ100Way) come with any
  account; Professional-tier (PrimateAI-3D, SpliceAI, COSMIC, OMIM) need the
  matching license.
- Data versions after the 2026-07-09 refresh: Ensembl 115, RefSeq
  GCF_000001405.40-RS_2025_08, ClinVar 20260328, dbSNP 157, COSMIC 103,
  gnomAD 4.1, plus added DANN / Gerp / AlphaMissense / PromoterAI.
