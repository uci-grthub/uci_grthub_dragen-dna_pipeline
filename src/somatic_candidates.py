#!/usr/bin/env python3
"""Genome-wide candidate variant list for the tumor-only client deliverable.

Less curated than panel_client_summary.py: every rare protein-altering PASS call
that DRAGEN did not tag as population germline, across the whole genome, with
panel genes flagged. Adds a sheet with every panel-gene call (any impact or
FILTER) and why it was left out of the summary, plus a filter-funnel sheet.

Usage: somatic_candidates.py CSQ_TSV ANNOT_TSV COUNTS_TSV PANEL_DB OUT_XLSX OUT_TSV PARAMS_YAML...
  CSQ_TSV    bcftools +split-vep -d rows (see rule candidates_table)
  ANNOT_TSV  nirvana_json_to_tsv.sh output for the merged candidate sites
  COUNTS_TSV sample <TAB> number of PASS non-Germline_DB records (step 1)
"""

import csv
import sqlite3
import sys

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill

from panel_client_summary import interpret, protein, read_params, short_clinvar

CSQ_COLS = [
    "chrom", "pos", "ref", "alt", "sample", "gt", "dp", "gene", "feature",
    "canonical", "consequence", "impact", "hgvsc", "hgvsp", "variant_id",
    "filter", "vaf", "germline_status", "cytoband",
]
ANNOT_COLS = [
    "variant_id", "rsid", "gnomad_af", "gnomad_hom", "clinvar", "clinvar_review",
    "clinvar_conditions", "revel", "phylop",
]
AF_MAX = 0.001
COMMON_AF = 0.01  # panel sheet: drop common polymorphisms
HEADER_FILL = PatternFill("solid", fgColor="DDDDDD")
PANEL_FILL = PatternFill("solid", fgColor="FFF2CC")


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def variant_class(gs):
    if gs and gs.startswith("Somatic"):
        return "somatic"
    if gs == "Germline_Proxi":
        return "unresolved"
    if gs == "Germline_DB":
        return "likely germline (population DB)"
    return "untagged"


def interpretation(r):
    if r["germline_status"] in ("", "."):
        return "untagged (MNV / complex record)"
    return interpret(r)


def clinvar(sig):
    c = short_clinvar(sig)
    return "" if c == "-" else c


def is_rare(af):
    return af is None or af < AF_MAX


def read_tsv(path, cols):
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            yield dict(zip(cols, row))


def write_sheet(ws, header, rows, widths, panel_col=None):
    ws.append(header)
    for c in ws[1]:
        c.font = Font(bold=True)
        c.fill = HEADER_FILL
    for row in rows:
        ws.append(row)
        if panel_col is not None and row[panel_col]:
            for c in ws[ws.max_row]:
                c.fill = PANEL_FILL
    for i, w in enumerate(widths):
        ws.column_dimensions[chr(65 + i) if i < 26 else "A" + chr(65 + i - 26)].width = w
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions


def main():
    csq_tsv, annot_tsv, counts_tsv, db_path, out_xlsx, out_tsv, *param_paths = sys.argv[1:]
    params = read_params(param_paths)
    samples = sorted(params)

    annot = {a["variant_id"]: a for a in read_tsv(annot_tsv, ANNOT_COLS)}
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    panel_cat = {
        r["gene_symbol"]: r["category"]
        for r in con.execute("SELECT gene_symbol, category FROM panel_genes")
    }

    # Filter funnel: distinct variant ids per sample after each step.
    steps = ["PASS, not Germline_DB", "on a canonical Ensembl transcript",
             "high / moderate impact", f"gnomAD AF < {AF_MAX} or absent"]
    funnel = {s: [set() for _ in steps] for s in samples}

    cand = {}
    for r in read_tsv(csq_tsv, CSQ_COLS):
        s, vid = r["sample"], r["variant_id"]
        if r["gt"] in ("0/0", "./.", ".", "0|0"):
            continue
        if not (r["feature"].startswith("ENST") and r["canonical"] == "True"):
            continue
        funnel[s][1].add(vid)
        if r["impact"] not in ("high", "moderate"):
            continue
        funnel[s][2].add(vid)
        a = annot.get(vid, {})
        af = num(a.get("gnomad_af"))
        if not is_rare(af):
            continue
        funnel[s][3].add(vid)
        r["variant_class"] = variant_class(r["germline_status"])
        key = (s, vid, r["gene"])
        if key in cand:
            continue
        cand[key] = [
            s, r["gene"], panel_cat.get(r["gene"], ""), r["chrom"], int(r["pos"]),
            r["ref"], r["alt"], r["consequence"].replace("&", ", "), r["impact"],
            r["hgvsc"].split(":")[-1], protein(r["hgvsp"], ""),
            num(r["vaf"]), int(r["dp"]) if r["dp"].isdigit() else None,
            r["germline_status"], interpretation(r), af,
            clinvar(a.get("clinvar", "")),
            num(a.get("revel")), a.get("rsid", ""), r["cytoband"],
        ]

    cand_rows = sorted(
        cand.values(), key=lambda x: (x[2] == "", x[1], x[0], x[3], x[4])
    )
    cand_header = [
        "sample", "gene", "panel_category", "chrom", "pos", "ref", "alt",
        "consequence", "impact", "HGVSc", "protein", "VAF", "depth",
        "GermlineStatus", "interpretation", "gnomAD_AF", "ClinVar", "REVEL",
        "rsID", "cytoband",
    ]

    # Panel genes, every non-population-germline call on the canonical transcript.
    panel_rows = []
    q = """
        SELECT DISTINCT sample, gene_symbol, gene_category, chrom, pos, ref, alt,
               consequence, impact, hgvsc, hgvsp, vaf, dp, filter, germline_status,
               variant_class, gnomad_af, clinvar_significance, revel, rsid
          FROM panel_annotated
         WHERE feature LIKE 'ENST%' AND canonical = 'True'
           AND variant_class <> 'likely germline (population DB)'
         ORDER BY gene_category, gene_symbol, sample, pos
    """
    for r in con.execute(q):
        r = dict(r)
        af = num(r["gnomad_af"])
        if af is not None and af >= COMMON_AF:
            continue  # common polymorphism (DRAGEN leaves failed calls untagged)
        why = []
        if r["filter"] != "PASS":
            why.append(f"DRAGEN filter: {r['filter']}")
        if r["impact"] not in ("high", "moderate"):
            why.append(f"{r['impact']} impact")
        if not is_rare(af):
            why.append(f"gnomAD AF {af:.3g}")
        panel_rows.append([
            r["sample"], r["gene_symbol"], r["gene_category"], r["chrom"], r["pos"],
            r["ref"], r["alt"], r["consequence"].replace("&", ", "), r["impact"],
            (r["hgvsc"] or "").split(":")[-1], protein(r["hgvsp"] if r["hgvsp"] != "." else "", ""),
            num(r["vaf"]), r["dp"], r["filter"], r["germline_status"],
            interpretation(r) if r["filter"] == "PASS" else "failed DRAGEN filter",
            af, clinvar(r["clinvar_significance"] or ""),
            num(r["revel"]), "; ".join(why) or "in summary",
        ])
    panel_header = [
        "sample", "gene", "panel_category", "chrom", "pos", "ref", "alt",
        "consequence", "impact", "HGVSc", "protein", "VAF", "depth", "FILTER",
        "GermlineStatus", "interpretation", "gnomAD_AF", "ClinVar", "REVEL",
        "excluded_from_summary_because",
    ]

    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cand_header)
        for row in cand_rows:
            w.writerow(["" if v is None else v for v in row])

    wb = Workbook()
    ws = wb.active
    ws.title = "README"
    p = params[samples[0]]
    cov = ", ".join(f"{s} {num(params[s].get('median_cov')) or 0:.0f}x" for s in samples)
    call = ", ".join(f"{s} {params[s].get('callability_thresh', '?')}" for s in samples)
    n_by = {s: sum(1 for r in cand_rows if r[0] == s) for s in samples}
    n_panel = {s: sum(1 for r in cand_rows if r[0] == s and r[2]) for s in samples}
    lines = [
        ("Candidate variant list, PT2 and PT5 organoids (tumor-only WGS)", True),
        ("", False),
        ("What this is", True),
        ("A broader list than the one-page panel summary: every rare protein-altering "
         "small variant called in each tumor across the whole genome. Panel genes "
         "(PI3K, RAS/MAPK, upstream RTK, subtype context) are flagged in column "
         "panel_category and shaded yellow.", False),
        ("", False),
        ("Filters (sheet Candidates)", True),
        ("1. DRAGEN 4.4.7 tumor-only call, FILTER = PASS.", False),
        ("2. Not tagged Germline_DB by DRAGEN (gnomAD / 1000 Genomes population allele).", False),
        ("3. Consequence on the canonical Ensembl transcript (Nirvana / Illumina Connected "
         "Annotations).", False),
        ("4. High or moderate impact (missense, nonsense, frameshift, in-frame indel, "
         "splice donor/acceptor, start/stop loss).", False),
        (f"5. gnomAD allele frequency below {AF_MAX} or absent.", False),
        (f"Rows: {', '.join(f'{s} {n_by[s]} ({n_panel[s]} in panel genes)' for s in samples)}. "
         "Sheet Counts shows how many variants remain after each step.", False),
        ("", False),
        ("Sheet Panel genes, all calls", True),
        ("Every call in the panel genes that is not a population germline allele, "
         "whatever its impact or DRAGEN filter. Column excluded_from_summary_because "
         "says why a call is not in the one-page summary. Common polymorphisms "
         f"(gnomAD AF >= {COMMON_AF}) are left out; DRAGEN does not germline-tag calls "
         "that fail its filters, so without this step they would fill the sheet.", False),
        ("", False),
        ("Settings", True),
        (f"Median autosomal coverage: {cov}. Callability threshold (reads): {call}. "
         f"Calibration: {p.get('calibration', '?')}.", False),
        ("", False),
        ("Limitations", True),
        ("No matched normal tissue: somatic status is inferred. GermlineStatus "
         "Germline_Proxi means DRAGEN judged the call germline from its VAF (~0.5); in "
         "organoids that are close to 100% tumor, clonal somatic variants also sit near "
         "0.5, so these calls are unresolved, not proven germline. Most rows in this "
         "list are Germline_Proxi; many are likely rare private germline variants.", False),
        ("Coverage is about 25x: variants below ~15% VAF have few supporting reads and "
         "may be missed.", False),
        ("Calls have not been individually reviewed in IGV.", False),
    ]
    ws.column_dimensions["A"].width = 110
    for text, bold in lines:
        ws.append([text])
        c = ws[ws.max_row][0]
        c.alignment = Alignment(wrap_text=True, vertical="top")
        if bold:
            c.font = Font(bold=True)

    widths = [7, 11, 14, 7, 11, 8, 8, 24, 9, 22, 22, 6, 6, 16, 34, 10, 22, 7, 14, 9]
    write_sheet(wb.create_sheet("Candidates"), cand_header, cand_rows, widths, panel_col=2)
    widths = [7, 10, 14, 7, 11, 8, 8, 24, 9, 22, 22, 6, 6, 22, 16, 30, 10, 22, 7, 40]
    write_sheet(wb.create_sheet("Panel genes, all calls"), panel_header, panel_rows, widths)

    ws = wb.create_sheet("Counts")
    n_raw = {r["sample"]: int(r["n"]) for r in read_tsv(counts_tsv, ["sample", "n"])}
    count_rows = [[steps[0]] + [n_raw[s] for s in samples]]
    count_rows += [[st] + [len(funnel[s][i]) for s in samples] for i, st in enumerate(steps) if i]
    write_sheet(ws, ["filter step (cumulative)"] + samples, count_rows, [40] + [10] * len(samples))

    wb.save(out_xlsx)
    for s in samples:
        print(s, n_raw[s], [len(x) for x in funnel[s][1:]], "rows", n_by[s], "panel", n_panel[s])


if __name__ == "__main__":
    main()
