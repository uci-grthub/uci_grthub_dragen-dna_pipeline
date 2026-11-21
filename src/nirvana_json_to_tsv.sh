#!/usr/bin/env bash
# Nirvana JSON (.json.gz) -> one TSV row per variant (vid):
# vid, dbSNP, gnomAD AF, gnomAD hom count, ClinVar classification, review status,
# conditions, REVEL, phyloP. Usage: nirvana_json_to_tsv.sh IN.json.gz > OUT.tsv
set -euo pipefail
zcat "$1" | jq -r '
    .positions[].variants[]
    | [ .vid,
        (.dbsnp // [] | join(",")),
        (.gnomad.allAf // ""),
        (.gnomad.allHc // ""),
        ([.clinvar[]?.classifications.germlineClassification.classification]
            | map(select(.)) | unique | join(";")),
        ([.clinvar[]?.classifications.germlineClassification.reviewStatus]
            | map(select(.)) | unique | join(";")),
        ([.clinvar[]?.classifications.germlineClassification.conditions[]?.traits[]?.name.value]
            | map(select(. and . != "not provided")) | unique | join(";")),
        (.revel.score // ""),
        (.phylopScore // "")
      ] | @tsv' \
    | awk -F'\t' '!seen[$1]++'
