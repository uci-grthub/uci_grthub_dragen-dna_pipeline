#!/usr/bin/env bash
set -euo pipefail

# Generate a sample CSV (sample,R1,R2) for dragen-run-exome.sh by pairing
# R1/R2 FASTQ files found in an input directory.
#
# Sample name is derived from the R1 filename by stripping the read/suffix,
# e.g.  xR089-L4-G4-P159-TCTACGCC-GTACGCCG-R1.fastq.gz
#   ->  sample = xR089-L4-G4-P159-TCTACGCC-GTACGCCG

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --input-dir PATH   Directory containing paired FASTQ files (default: data/FASTQ)
  --output PATH      Output CSV path (default: stdout)
  --r1-tag STR       R1 marker in filenames (default: R1)
  --r2-tag STR       R2 marker in filenames (default: R2)
  --abs              Write absolute paths instead of basenames
  -h, --help         Show this help
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="$SCRIPT_DIR/../data/FASTQ"
OUTPUT=""
R1_TAG="R1"
R2_TAG="R2"
ABS=false

while [ $# -gt 0 ]; do
    case "$1" in
        --input-dir) INPUT_DIR="$2"; shift 2;;
        --output) OUTPUT="$2"; shift 2;;
        --r1-tag) R1_TAG="$2"; shift 2;;
        --r2-tag) R2_TAG="$2"; shift 2;;
        --abs) ABS=true; shift;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown option: $1" >&2; usage; exit 1;;
    esac
done

if [ ! -d "$INPUT_DIR" ]; then
    echo "Input dir not found: $INPUT_DIR" >&2
    exit 1
fi

INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"

emit() {
    printf '# sample,R1,R2\n'
    local found=0
    shopt -s nullglob
    for r1 in "$INPUT_DIR"/*"${R1_TAG}"*.fastq.gz "$INPUT_DIR"/*"${R1_TAG}"*.fq.gz; do
        local base r2 sample
        base="$(basename "$r1")"
        r2="${r1/${R1_TAG}/${R2_TAG}}"
        if [ ! -e "$r2" ]; then
            echo "No R2 mate for $base (expected $(basename "$r2")); skipping" >&2
            continue
        fi
        # sample = filename up to and including the token before the R1 tag
        sample="${base%%-${R1_TAG}.*}"
        sample="${sample%%_${R1_TAG}_*}"
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
