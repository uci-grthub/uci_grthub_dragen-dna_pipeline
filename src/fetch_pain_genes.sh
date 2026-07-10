#!/usr/bin/env bash
# Build a sourced, reproducible pain-gene panel (data/pain_genes.tsv).
#
# Two provenance sources, tagged per gene in the `source` column:
#   GO   - pain biology: genes annotated (QuickGO, human Swiss-Prot) to
#          "sensory perception of pain" (GO:0019233), "response to pain"
#          (GO:0048265), "behavioral response to pain" (GO:0048266), and
#          "opioid receptor signaling pathway" (GO:0038003), incl. descendants.
#   CPIC / PharmGKB - pain management pharmacogenomics: analgesic gene-drug
#          pairs from CPIC clinical guidelines (+ a few PharmGKB VIP genes).
#
# Ensembl gene ids are joined from the configured GENCODE GTF.
# Requires network (QuickGO). Re-run to refresh; the GO half is versioned by
# download date (recorded in the file header comment).
#
# Usage: fetch_pain_genes.sh [gene_gtf.gz] [out_tsv]
#   gene_gtf.gz defaults to config/config.yaml:gene_gtf; out_tsv to data/pain_genes.tsv
set -euo pipefail

GTF=${1:-$(grep -E '^gene_gtf:' config/config.yaml | sed 's/[^"]*"//;s/".*//')}
OUT=${2:-data/pain_genes.tsv}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

GO_TERMS="GO:0019233,GO:0048265,GO:0048266,GO:0038003"
QUICKGO="https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"

echo "[fetch] QuickGO pain-biology annotations ($GO_TERMS)" >&2
curl -sSL --max-time 120 -H 'Accept:text/tsv' \
  "${QUICKGO}?goId=${GO_TERMS}&goUsage=descendants&taxonId=9606&geneProductType=protein&geneProductSubset=Swiss-Prot" \
  -o "$WORK/go.tsv"

# symbol -> semicolon list of distinct GO term ids (evidence)
awk -F'\t' 'NR>1 && $3!="" {
    key=$3; go=$5
    if (!(key SUBSEP go in seen)) { seen[key SUBSEP go]=1;
        ev[key] = (key in ev ? ev[key]";"go : go) }
} END { for (g in ev) print g"\tGO\t"ev[g] }' "$WORK/go.tsv" \
  | sort > "$WORK/go_genes.tsv"

echo "[fetch] CPIC / PharmGKB analgesic pharmacogenes (pinned)" >&2
# gene <tab> source <tab> evidence(guideline: drugs)
cat > "$WORK/pgx.tsv" <<'PGX'
CYP2D6	CPIC	CPIC opioids: codeine,tramadol,hydrocodone,oxycodone
OPRM1	CPIC	CPIC opioids: response
COMT	CPIC	CPIC opioids: response
CYP2C9	CPIC	CPIC NSAIDs: celecoxib,flurbiprofen,ibuprofen,meloxicam,piroxicam
UGT2B7	PharmGKB	PharmGKB: morphine glucuronidation
ABCB1	PharmGKB	PharmGKB: opioid transport (P-glycoprotein)
CYP2B6	PharmGKB	PharmGKB: methadone metabolism
PGX
sort -o "$WORK/pgx.tsv" "$WORK/pgx.tsv"

echo "[merge] union GO + PGx, tag combined source" >&2
# join on gene symbol (full outer): combine source + evidence
join -t $'\t' -a1 -a2 -e '' -o '0,1.2,1.3,2.2,2.3' \
     -j1 "$WORK/go_genes.tsv" "$WORK/pgx.tsv" 2>/dev/null \
  | awk -F'\t' '{
      src = ($2!="" && $4!="") ? $2";"$4 : ($2!="" ? $2 : $4)
      ev  = ($3!="" && $5!="") ? $3" | "$5 : ($3!="" ? $3 : $5)
      print $1"\t"src"\t"ev
    }' | sort > "$WORK/merged.tsv"

echo "[join] Ensembl gene ids from GTF" >&2
# symbol -> ENSG (strip version) from GENCODE gene lines
zcat "$GTF" | gawk -F'\t' '$3=="gene" {
    name=$9; sub(/.*gene_name "/,"",name); sub(/".*/,"",name)
    id=$9;   sub(/.*gene_id "/,"",id);     sub(/".*/,"",id); sub(/\..*/,"",id)
    if (!(name in G)) G[name]=id
    print name"\t"id
}' | sort -u > "$WORK/ensg.tsv"

# emit final TSV
DATE=$(date +%Y-%m-%d)
{
  echo "# pain-gene panel — generated $DATE by src/fetch_pain_genes.sh"
  echo "# sources: GO (QuickGO ${GO_TERMS}, human Swiss-Prot, +descendants); CPIC/PharmGKB analgesic pharmacogenes"
  echo -e "gene_symbol\tsource\tevidence\tensembl_gene_id"
  join -t $'\t' -a1 -e '' -o '0,1.2,1.3,2.2' -j1 "$WORK/merged.tsv" \
       <(awk -F'\t' '!s[$1]++' "$WORK/ensg.tsv") 2>/dev/null \
    | sort -k1,1
} > "$OUT"

N=$(grep -vc '^#\|^gene_symbol' "$OUT")
echo "[done] wrote $OUT — $N genes" >&2
