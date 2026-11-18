#!/usr/bin/env bash
set -euo pipefail

# Generate a sample CSV (sample,R1,R2) for the Snakemake workflow by pairing
# R1/R2 FASTQ files found in an input directory.
#
# Default sample name is derived from the R1 filename by stripping the
# read/suffix, e.g.  xR089-L4-G4-P159-TCTACGCC-GTACGCCG-R1.fastq.gz
#   ->  sample = xR089-L4-G4-P159-TCTACGCC-GTACGCCG
#
# With --metadata, the sample name is instead looked up in a metadata CSV
# (columns: Sample,Index,i7 Barcode,i5 Barcode (MiSeq),i5 Barcode (NovaSeq v1.5))
# by matching "-<i7>-<i5 NovaSeq>-" in the FASTQ filename. FASTQs with no
# matching barcode pair are skipped with a warning.

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --input-dir PATH   Directory containing paired FASTQ files (default: data/FASTQ)
  --output PATH      Output CSV path (default: stdout)
  --r1-tag STR       R1 marker in filenames (default: R1)
  --r2-tag STR       R2 marker in filenames (default: R2)
  --abs              Write absolute paths instead of basenames
  --metadata PATH    Metadata CSV mapping i7/i5 barcodes to sample names
  -h, --help         Show this help
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="$SCRIPT_DIR/../data/FASTQ"
OUTPUT=""
R1_TAG="R1"
R2_TAG="R2"
ABS=false
METADATA=""

while [ $# -gt 0 ]; do
    case "$1" in
        --input-dir) INPUT_DIR="$2"; shift 2;;
        --output) OUTPUT="$2"; shift 2;;
        --r1-tag) R1_TAG="$2"; shift 2;;
        --r2-tag) R2_TAG="$2"; shift 2;;
        --abs) ABS=true; shift;;
        --metadata) METADATA="$2"; shift 2;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown option: $1" >&2; usage; exit 1;;
    esac
done

if [ ! -d "$INPUT_DIR" ]; then
    echo "Input dir not found: $INPUT_DIR" >&2
    exit 1
fi

INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"

# barcode key "<i7>-<i5 NovaSeq v1.5>" -> sample name
declare -A BARCODE_SAMPLE=()
if [ -n "$METADATA" ]; then
    if [ ! -f "$METADATA" ]; then
        echo "Metadata file not found: $METADATA" >&2
        exit 1
    fi
    while IFS=, read -r name _idx i7 _i5miseq i5nova _rest; do
        name="${name//[[:space:]]/}"; i7="${i7//[[:space:]]/}"
        i5nova="${i5nova//[[:space:]]/}"
        [ -z "$name" ] || [ "$name" = "Sample" ] && continue
        [ -z "$i7" ] || [ -z "$i5nova" ] && continue
        BARCODE_SAMPLE["$i7-$i5nova"]="$name"
    done < <(tr -d '\r' < "$METADATA")
    if [ "${#BARCODE_SAMPLE[@]}" -eq 0 ]; then
        echo "No samples parsed from metadata: $METADATA" >&2
        exit 1
    fi
fi

emit() {
    printf '# sample,R1,R2\n'
    local found=0
    shopt -s nullglob
    # read tag must be delimited ("-R1." or "_R1_"), so run IDs such as
    # "xR128" are not mistaken for the read marker
    for r1 in "$INPUT_DIR"/*-"${R1_TAG}".f*q.gz "$INPUT_DIR"/*_"${R1_TAG}"_*.f*q.gz; do
        local base r2 sample
        base="$(basename "$r1")"
        if [[ "$base" == *-"${R1_TAG}".* ]]; then
            r2="$INPUT_DIR/${base/-${R1_TAG}./-${R2_TAG}.}"
        else
            r2="$INPUT_DIR/${base/_${R1_TAG}_/_${R2_TAG}_}"
        fi
        if [ ! -e "$r2" ]; then
            echo "No R2 mate for $base (expected $(basename "$r2")); skipping" >&2
            continue
        fi
        # sample = filename up to and including the token before the R1 tag
        sample="${base%%-${R1_TAG}.*}"
        sample="${sample%%_${R1_TAG}_*}"
        if [ -n "$METADATA" ]; then
            local key matched=""
            for key in "${!BARCODE_SAMPLE[@]}"; do
                if [[ "$base" == *"-${key}-"* ]]; then
                    matched="${BARCODE_SAMPLE[$key]}"
                    break
                fi
            done
            if [ -z "$matched" ]; then
                echo "No metadata barcode match for $base; skipping" >&2
                continue
            fi
            sample="$matched"
        fi
        if [ "$ABS" = true ]; then
            printf '%s,%s,%s\n' "$sample" "$r1" "$r2"
        else
            printf '%s,%s,%s\n' "$sample" "$(basename "$r1")" "$(basename "$r2")"
        fi
        found=$((found+1))
    done
    if [ "$found" -eq 0 ]; then
        echo "No R1 FASTQ files matched in $INPUT_DIR" >&2
        return 1
    fi
}

if [ -n "$OUTPUT" ]; then
    emit > "$OUTPUT"
    echo "Wrote sample file: $OUTPUT" >&2
else
    emit
fi
