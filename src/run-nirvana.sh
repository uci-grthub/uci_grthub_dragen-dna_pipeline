#!/usr/bin/env bash
set -euo pipefail

# Standalone Nirvana annotation pass over DRAGEN output already on disk.
# Use when a run finished with variant calling done but Nirvana skipped
# (e.g. a broken NIRVANA_BIN path). Annotates every hard-filtered VCF found
# under <output-base>/<sample>/<sample>.hard-filtered.vcf.gz.

usage() {
    cat <<EOF
Usage: $(basename "$0") --output-base PATH [options]

Options:
  --output-base PATH   Base output dir containing per-sample subdirectories (required)
  --vcf-suffix STR     VCF suffix to look for (default: .hard-filtered.vcf.gz)
  --nirvana-bin PATH   Nirvana binary (default: DRAGEN 4.4.7 bundled)
  --cache PATH         Nirvana cache prefix (default: GRCh38 Both)
  --dat PATH           Nirvana reference .dat (default: GRCh38)
  --sd PATH            Supplementary annotation dir (default: GRCh38)
  --force              Re-annotate even if a .json.gz output already exists
  --dry-run            Print commands but do not execute
  -h, --help           Show this help
EOF
}

OUTPUT_BASE=""
VCF_SUFFIX=".hard-filtered.vcf.gz"
NIRVANA_BIN="/opt/dragen/4.4.7/share/nirvana/Nirvana"
NIRVANA_CACHE="/staging/human/reference/hg38/annotation_data_files/Cache/GRCh38"
NIRVANA_DAT="/staging/human/reference/hg38/annotation_data_files/References/Homo_sapiens.GRCh38.Nirvana.dat"
NIRVANA_SD="/staging/human/reference/hg38/annotation_data_files/SupplementaryAnnotation/GRCh38_curated"
FORCE=false
DRY_RUN=false

if [ $# -eq 0 ]; then usage; exit 1; fi

while [ $# -gt 0 ]; do
    case "$1" in
        --output-base) OUTPUT_BASE="$2"; shift 2;;
        --vcf-suffix) VCF_SUFFIX="$2"; shift 2;;
        --nirvana-bin) NIRVANA_BIN="$2"; shift 2;;
        --cache) NIRVANA_CACHE="$2"; shift 2;;
        --dat) NIRVANA_DAT="$2"; shift 2;;
        --sd) NIRVANA_SD="$2"; shift 2;;
        --force) FORCE=true; shift;;
        --dry-run) DRY_RUN=true; shift;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown option: $1" >&2; usage; exit 1;;
    esac
done

if [ -z "$OUTPUT_BASE" ]; then
    echo "--output-base is required" >&2
    exit 1
fi
if [ ! -x "$NIRVANA_BIN" ]; then
    echo "Nirvana binary not found or not executable: $NIRVANA_BIN" >&2
    exit 1
fi

shopt -s nullglob
found=0
for vcf in "$OUTPUT_BASE"/*/*"$VCF_SUFFIX"; do
    found=$((found+1))
    dir="$(dirname "$vcf")"
    base="$(basename "$vcf" "$VCF_SUFFIX")"
    out="$dir/$base"

    if [ "$FORCE" = false ] && [ -e "$out.json.gz" ]; then
        echo "Already annotated, skipping: $out.json.gz (use --force to redo)" >&2
        continue
    fi

    cmd=("$NIRVANA_BIN" -c "$NIRVANA_CACHE" -r "$NIRVANA_DAT" --sd "$NIRVANA_SD" -i "$vcf" -o "$out")
    if [ "$DRY_RUN" = true ]; then
        printf "DRY: %s\n" "${cmd[*]}"
    else
        echo "Annotating: $vcf"
        "${cmd[@]}"
    fi
done

if [ "$found" -eq 0 ]; then
    echo "No VCFs matching *$VCF_SUFFIX under $OUTPUT_BASE" >&2
    exit 1
fi

echo "Done."
