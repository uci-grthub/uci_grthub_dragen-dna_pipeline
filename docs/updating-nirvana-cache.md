# Updating the Nirvana / ICA annotation cache

The bundled Illumina Connected Annotations (ICA) binary is **3.26.0** (under
`/opt/dragen/4.4.7/share/nirvana/`). The transcript cache currently on disk
(`.../Cache/GRCh38/Both.transcripts.ndb`) is an **old-format cache that ICA 3.26
silently rejects** — it "loads" instantly but contributes no transcript data
source, so annotated output has **no gene / HGVS / consequence** fields
(gnomAD, ClinVar, dbSNP, SpliceAI, REVEL etc. from the SD dir still work).

Fixing this means downloading an ICA-3.26-compatible cache with `DataManager`.

## Why it is currently blocked

`DataManager make-config` and `download` both require Illumina credentials at:

```
~/.ilmnAnnotations/credentials.json
```

That file does not exist on this machine, so the download cannot proceed until
you provide it. It is tied to your Illumina Connected Annotations account /
license (contact `annotation_support@illumina.com` if you don't have one).

## How to tell the cache is broken

Run ICA on any VCF and look at the **Data Source** table it prints:

```bash
N=/opt/dragen/4.4.7/share/nirvana/Nirvana
$N -c <cache-dir> -r <ref.dat> --sd <sd-dir> -i some.vcf.gz -o /tmp/probe
```

- **Broken:** the table lists only `Gene`-type rows (clingen, cosmic, gnomad,
  omim) and **no `Transcript` row**. `grep -c '"transcripts"'` on the JSON = 0.
- **Working:** you also see rows like `RefSeq  Transcript  Loaded` /
  `Ensembl  Transcript  Loaded`, and variants carry a `transcripts` array.

## Steps to update

1. **Install credentials.** Obtain `credentials.json` from your Illumina ICA
   account and place it here:

   ```bash
   mkdir -p ~/.ilmnAnnotations
   cp /path/to/credentials.json ~/.ilmnAnnotations/credentials.json
   chmod 600 ~/.ilmnAnnotations/credentials.json
   ```

2. **Generate the download config** for the assembly:

   ```bash
   DM=/opt/dragen/4.4.7/share/nirvana/DataManager
   $DM make-config --ref GRCh38 -o ~/.ilmnAnnotations
   # writes ~/.ilmnAnnotations/GRCh38_annotation_config.json
   ```

3. **Download** into a fresh directory (large; needs network):

   ```bash
   NEW=/staging/human/reference/hg38/annotation_data_files_ica326
   mkdir -p "$NEW"
   $DM download --ref GRCh38 --dir "$NEW" --thread 4
   ```

   This populates a matching set:
   - `Cache/GRCh38/`               transcript cache (the missing piece)
   - `References/…Nirvana.dat`     reference sequence
   - `SupplementaryAnnotation/GRCh38/`  gnomAD, ClinVar, dbSNP, SpliceAI, …

4. **Verify** the new cache loads a transcript source:

   ```bash
   $N -c "$NEW/Cache/GRCh38" \
      -r "$NEW/References/Homo_sapiens.GRCh38.Nirvana.dat" \
      --sd "$NEW/SupplementaryAnnotation/GRCh38" \
      -i output/<sample>/<sample>.hard-filtered.vcf.gz -o /tmp/probe
   zcat /tmp/probe.json.gz | grep -c '"transcripts"'   # should be > 0
   ```

5. **Repoint the pipeline** in `config/config.yaml`:

   ```yaml
   nirvana_cache: "/staging/human/reference/hg38/annotation_data_files_ica326/Cache/GRCh38"
   nirvana_dat:   "/staging/human/reference/hg38/annotation_data_files_ica326/References/Homo_sapiens.GRCh38.Nirvana.dat"
   nirvana_sd:    "/staging/human/reference/hg38/annotation_data_files_ica326/SupplementaryAnnotation/GRCh38"
   ```

   Update the same three `NIRVANA_*` defaults in `src/run-nirvana.sh` and
   `src/dragen-run-exome.sh` if you use those directly.

6. **Re-annotate:**

   ```bash
   bash src/run-nirvana.sh --output-base output --force      # standalone
   # or
   pixi run run                                              # via Snakemake
   ```

## Notes

- The freshly downloaded `SupplementaryAnnotation/GRCh38` will have complete
  `.idx` sidecars, so the `GRCh38_curated` workaround (which drops DANN / GERP /
  PhyloP-primate) is no longer needed — point `nirvana_sd` at the real dir.
- `credentials.json` is a secret — keep it `chmod 600`, never commit it.
- Basic-tier sources (gnomAD, ClinVar, dbSNP, GERP, DANN, phyloP, TOPMed,
  1000 Genomes, DECIPHER, FusionCatcher, MITOMAP, MultiZ100Way) come with any
  account; Professional-tier (PrimateAI-3D, SpliceAI, COSMIC, OMIM) need the
  matching license.
