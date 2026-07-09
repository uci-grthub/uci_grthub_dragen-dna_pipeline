#!/usr/bin/env bash
set -euo pipefail

# Exome adaptation of dragen-run.sh.
# Runs DRAGEN over multiple exome samples, restricting variant calling to the
# IDT xGen Exome Research Panel v2 target BED (hg38) and reporting on-target
# coverage QC against the same panel.

VC_FILTER_DEFAULT=''
# VC_FILTER_DEFAULT='SNP filter:snp:QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0; INDEL filter:indel:QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'

# IDT xGen Exome Research Panel v2 BEDs (hg38)
EXOME_BED_DIR="/staging/human/reference/hg38/exome_beds/idt_x_gen_exome_research_panel_v2"
TARGET_BED_DEFAULT="$EXOME_BED_DIR/xgen-exome-hyb-panel-v2-targets-hg38.bed"
PROBE_BED_DEFAULT="$EXOME_BED_DIR/xgen-exome-hyb-panel-v2-probes-hg38.bed"
# Spike-in panel BEDs (extra capture regions added on top of the base kit)
SPIKE_TARGET_BED_DEFAULT="$EXOME_BED_DIR/spike-in-xgen-exome-v2-targets-hg38.bed"
SPIKE_PROBE_BED_DEFAULT="$EXOME_BED_DIR/spike-in-xgen-exome-v2-probes-hg38.bed"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Generalized driver for running DRAGEN on multiple exome samples.

Options:
  --reference PATH        DRAGEN reference dir (default: hg38 pangenome, hashtable v11)
  --input-dir PATH        Directory where FASTQ files live (optional when sample file contains full paths)
  --sample-file PATH      CSV file with lines: sample,R1,R2  (R1/R2 may be basenames or full paths)
  --output-base PATH      Base output directory (required)
  --target-bed PATH       Exome target BED for variant calling (default: xGen v2 targets hg38)
  --target-padding N      Bases to pad each side of target regions for variant calling (default: 100)
  --probe-bed PATH        Exome probe/capture BED for coverage QC (default: xGen v2 probes hg38)
  --with-spike-in         Merge the spike-in panel BEDs into the kit target/probe BEDs
                          (union). Use when the library prep included the spike-in.
  --spike-target-bed PATH Spike-in target BED (default: spike-in xGen v2 targets hg38)
  --spike-probe-bed PATH  Spike-in probe BED (default: spike-in xGen v2 probes hg38)
  --rgid STR              Read group ID (default: Illumina_RGID)
  --vc-filter STR         Variant-caller hard filter string
  --enable-variant-caller true|false  (default: true)
  --dry-run               Print commands but do not execute
  --create-only           Only create output directories and print commands (do not run DRAGEN)
  -h, --help              Show this help

Examples:
  $(basename "$0") --input-dir /data/FASTQ --sample-file data/samples.csv --output-base output

Sample file format:
  CSV with either 3 columns: sample,R1,R2  (R1/R2 can be basenames or full paths)
  Lines starting with # are ignored.

EOF
}

# defaults
INPUT_DIR=""
SAMPLE_FILE=""
REFERENCE="/staging/human/reference/hg38/pangenome_reference"
OUTPUT_BASE=""
TARGET_BED="$TARGET_BED_DEFAULT"
TARGET_PADDING=100
PROBE_BED="$PROBE_BED_DEFAULT"
WITH_SPIKE_IN=false
SPIKE_TARGET_BED="$SPIKE_TARGET_BED_DEFAULT"
SPIKE_PROBE_BED="$SPIKE_PROBE_BED_DEFAULT"
RGID="Illumina_RGID"
VC_FILTER="$VC_FILTER_DEFAULT"
ENABLE_VC=true
DRY_RUN=false
CREATE_ONLY=false
# Nirvana annotation defaults (hg38 / GRCh38; override via environment if needed)
NIRVANA_BIN="/opt/dragen/4.4.7/share/nirvana/Nirvana"
NIRVANA_CACHE="/staging/human/reference/hg38/annotation_data_files/Cache/GRCh38"
NIRVANA_DAT="/staging/human/reference/hg38/annotation_data_files/References/Homo_sapiens.GRCh38.Nirvana.dat"
NIRVANA_SD="/staging/human/reference/hg38/annotation_data_files/SupplementaryAnnotation/GRCh38_curated"

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --reference) REFERENCE="$2"; shift 2;;
        --input-dir) INPUT_DIR="$2"; shift 2;;
        --sample-file) SAMPLE_FILE="$2"; shift 2;;
        --output-base) OUTPUT_BASE="$2"; shift 2;;
        --target-bed) TARGET_BED="$2"; shift 2;;
        --target-padding) TARGET_PADDING="$2"; shift 2;;
        --probe-bed) PROBE_BED="$2"; shift 2;;
        --with-spike-in) WITH_SPIKE_IN=true; shift;;
        --spike-target-bed) SPIKE_TARGET_BED="$2"; shift 2;;
        --spike-probe-bed) SPIKE_PROBE_BED="$2"; shift 2;;
        --rgid) RGID="$2"; shift 2;;
        --vc-filter) VC_FILTER="$2"; shift 2;;
        --enable-variant-caller) ENABLE_VC="$2"; shift 2;;
        --dry-run) DRY_RUN=true; shift;;
        --create-only) CREATE_ONLY=true; shift;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown option: $1" >&2; usage; exit 1;;
    esac
done

if [ -z "$REFERENCE" ]; then
    echo "--reference is required" >&2
    exit 1
fi
if [ -z "$OUTPUT_BASE" ]; then
    echo "--output-base is required" >&2
    exit 1
fi
if [ -z "$SAMPLE_FILE" ]; then
    echo "--sample-file is required" >&2
    exit 1
fi
if [ ! -e "$TARGET_BED" ]; then
    echo "Target BED not found: $TARGET_BED" >&2
    exit 1
fi

# merge_beds A B OUT : union of two BEDs (sorted, overlaps/adjacent collapsed)
merge_beds() {
    local a="$1" b="$2" out="$3"
    if ! command -v bedtools >/dev/null 2>&1; then
        echo "bedtools not found on PATH; required for --with-spike-in" >&2
        exit 1
    fi
    cat "$a" "$b" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > "$out"
}

# When the prep used the spike-in panel, call variants + report coverage over
# the union of kit + spike-in regions.
if [ "$WITH_SPIKE_IN" = true ]; then
    for f in "$SPIKE_TARGET_BED" "$SPIKE_PROBE_BED"; do
        if [ ! -e "$f" ]; then
            echo "Spike-in BED not found: $f" >&2
            exit 1
        fi
    done
    mkdir -p "$OUTPUT_BASE"
    MERGED_TARGET="$OUTPUT_BASE/combined-targets-hg38.bed"
    MERGED_PROBE="$OUTPUT_BASE/combined-probes-hg38.bed"
    if [ "$DRY_RUN" = true ] || [ "$CREATE_ONLY" = true ]; then
        echo "DRY: merge $TARGET_BED + $SPIKE_TARGET_BED -> $MERGED_TARGET" >&2
        echo "DRY: merge $PROBE_BED + $SPIKE_PROBE_BED -> $MERGED_PROBE" >&2
    else
        echo "Merging kit + spike-in target BEDs -> $MERGED_TARGET" >&2
        merge_beds "$TARGET_BED" "$SPIKE_TARGET_BED" "$MERGED_TARGET"
        echo "Merging kit + spike-in probe BEDs -> $MERGED_PROBE" >&2
        merge_beds "$PROBE_BED" "$SPIKE_PROBE_BED" "$MERGED_PROBE"
    fi
    # use the combined BEDs for variant calling + coverage QC
    TARGET_BED="$MERGED_TARGET"
    PROBE_BED="$MERGED_PROBE"
fi

run_one() {
    local sample="$1"
    local r1="$2"
    local r2="$3"

    local outdir="$OUTPUT_BASE/$sample"
    mkdir -p "$outdir"

    local prefix="$sample"

    local cmd=(dragen -f \
        -r "$REFERENCE" \
        -1 "$r1" \
        -2 "$r2" \
        --enable-variant-caller "$ENABLE_VC" \
        --enable-map-align true \
        --enable-map-align-output true \
        --output-format BAM \
        --vc-target-bed "$TARGET_BED" \
        --vc-target-bed-padding "$TARGET_PADDING" \
        --RGID "$RGID" \
        --RGSM "$sample" \
        --output-directory "$outdir" \
        --output-file-prefix "$prefix" \
	--remove-duplicates true \
        --read-trimmers adapter,quality \
	--trim-adapter-read1 /staging/human/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa \
	--trim-adapter-read2 /staging/human/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa \
	--trim-min-quality 15 \
	--trim-min-length 20 \
        --qc-coverage-region-1 "$TARGET_BED" \
        --qc-coverage-reports-1 cov_report)

    # apply a hard filter only when one is provided (default empty; set via --vc-filter)
    if [ -n "$VC_FILTER" ]; then
        cmd+=(--vc-hard-filter "$VC_FILTER")
    fi

    # add probe/capture BED as a second coverage region when present
    if [ -e "$PROBE_BED" ]; then
        cmd+=(--qc-coverage-region-2 "$PROBE_BED" --qc-coverage-reports-2 cov_report)
    fi

    if [ "$DRY_RUN" = true ] || [ "$CREATE_ONLY" = true ]; then
        printf "DRY: %s\n" "${cmd[*]}"
    else
        echo "Running: ${cmd[*]}"
        "${cmd[@]}"
    fi

    # Attempt to annotate the hard-filtered VCF with Nirvana (if present)
    local vcf="$outdir/${prefix}.hard-filtered.vcf.gz"

    # construct Nirvana command using defaults (these can be overridden via env)
    local nirvana_cmd=("$NIRVANA_BIN" -c "$NIRVANA_CACHE" -r "$NIRVANA_DAT" --sd "$NIRVANA_SD" -i "$vcf" -o "$outdir/${prefix}")

    if [ -x "$NIRVANA_BIN" ] && [ -e "$vcf" ]; then
        if [ "$DRY_RUN" = true ] || [ "$CREATE_ONLY" = true ]; then
            printf "DRY: %s\n" "${nirvana_cmd[*]}"
        else
            echo "Running: ${nirvana_cmd[*]}"
            "${nirvana_cmd[@]}"
        fi
    else
        if [ ! -x "$NIRVANA_BIN" ]; then
            echo "Nirvana binary not found or not executable at $NIRVANA_BIN; skipping Nirvana for $sample" >&2
        elif [ ! -e "$vcf" ]; then
            echo "VCF not found: $vcf; skipping Nirvana for $sample" >&2
        fi
    fi
}

find_fastq_pair() {
    local sample="$1"
    # try to find R1 and R2 by globbing in INPUT_DIR
    local r1match=("$INPUT_DIR/${sample}"*R1*fastq* "$INPUT_DIR/${sample}"*R1*"")
    local r2match=("$INPUT_DIR/${sample}"*R2*fastq* "$INPUT_DIR/${sample}"*R2*"")
    # pick first existing
    local r1=""
    local r2=""
    for f in "${r1match[@]}"; do
        [ -e "$f" ] && { r1="$f"; break; }
    done
    for f in "${r2match[@]}"; do
        [ -e "$f" ] && { r2="$f"; break; }
    done
    printf '%s\t%s' "$r1" "$r2"
}

while IFS= read -r line || [ -n "$line" ]; do
    # skip comments and blank
    line_trimmed="${line%%#*}"
    line_trimmed="$(echo "$line_trimmed" | tr -d '\r' | sed -e 's/^\s*//' -e 's/\s*$//')"
    [ -z "$line_trimmed" ] && continue

    IFS=',' read -r sample col2 col3 <<< "$line_trimmed"
    sample="$(echo "$sample" | tr -d '"' | sed -e 's/^\s*//' -e 's/\s*$//')"

    r1=""
    r2=""
    if [ -n "${col2:-}" ] && [ -n "${col3:-}" ]; then
        r1="$col2"
        r2="$col3"
        # if not absolute, prefix with INPUT_DIR when provided
        if [[ "$r1" != /* ]] && [ -n "$INPUT_DIR" ]; then
            r1="$INPUT_DIR/$r1"
        fi
        if [[ "$r2" != /* ]] && [ -n "$INPUT_DIR" ]; then
            r2="$INPUT_DIR/$r2"
        fi
    else
        if [ -z "$INPUT_DIR" ]; then
            echo "Cannot infer R1/R2 for sample '$sample' without --input-dir or columns in sample file" >&2
            exit 1
        fi
        # attempt to find
        read -r r1 r2 < <(find_fastq_pair "$sample")
        if [ -z "$r1" ] || [ -z "$r2" ]; then
            echo "Fastq pair not found for sample $sample in $INPUT_DIR" >&2
            exit 1
        fi
    fi

    run_one "$sample" "$r1" "$r2"
done < "$SAMPLE_FILE"

echo "Done."
