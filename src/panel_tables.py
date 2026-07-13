#!/usr/bin/env python3
"""Flatten the curated panel CSV into two TSVs for loading into the panel DB.

  genes.tsv    - one row per panel gene, all CSV columns verbatim.
  keyvars.tsv  - one row per (gene, rsID) parsed out of the `key_variants` column,
                 so a-priori variants can be flagged in the cohort by exact rsID.

Only rsIDs can be matched exactly. Genes whose key variants are described as star
alleles, haplotypes or a variant class ("LoF variants", "GWAS locus") yield no rows
here and are therefore never flagged — that is a limitation of the panel text, not
evidence the variant is absent.

Usage: panel_tables.py panel.csv genes.tsv keyvars.tsv
"""
import csv
import re
import sys

csv_path, genes_path, keyvars_path = sys.argv[1], sys.argv[2], sys.argv[3]

GENE_COLS = [
    "gene_symbol",
    "ensembl_gene_id",
    "category",
    "priority",
    "mechanism",
    "key_variants",
    "evidence_note",
    "exome_capture_note",
]


def clean(s):
    # Tabs/newlines would break the sqlite .import; the CSV prose has neither,
    # but a stray one would silently shift every column.
    return " ".join((s or "").split())


with open(csv_path, newline="") as fh:
    rows = [r for r in csv.DictReader(fh) if (r.get("gene_symbol") or "").strip()]

with open(genes_path, "w", newline="") as gf, open(keyvars_path, "w", newline="") as kf:
    genes = csv.writer(gf, delimiter="\t", lineterminator="\n")
    keyvars = csv.writer(kf, delimiter="\t", lineterminator="\n")
    for r in rows:
        genes.writerow([clean(r.get(c)) for c in GENE_COLS])
        for rsid in re.findall(r"rs\d+", r.get("key_variants") or ""):
            keyvars.writerow([r["gene_symbol"], rsid, clean(r.get("key_variants"))])
