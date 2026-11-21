#!/usr/bin/env python3
"""Short client summary of a tumor-only gene-panel DB.

Writes a 3-sheet xlsx (gene matrix, alterations, methods) and a one-page
Markdown summary from the panel SQLite DB built by rule panel_db.

usage: panel_client_summary.py DB OUT_XLSX OUT_MD PARAMS_YAML [PARAMS_YAML ...]

PARAMS_YAML are the per-sample <sample>.somatic_params.yaml files of the
calibration the panel reads (median coverage, callability threshold, CNV bin
width go into the methods text).
"""

import re
import sqlite3
import sys

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill

# Same criteria as the pathway_hits view, without the FILTER = PASS clause,
# so DRAGEN-filtered calls can be listed separately.
CANDIDATES = """
SELECT gene_symbol, gene_category, sample, chrom, pos, ref, alt, dp, vaf,
       filter, germline_status, variant_class, consequence, hgvsc, hgvsp,
       gnomad_af, clinvar_significance
FROM panel_annotated
WHERE impact IN ('high', 'moderate')
  AND feature LIKE 'ENST%' AND canonical = 'True'
  AND variant_class <> 'likely germline (population DB)'
  AND (gnomad_af IS NULL OR gnomad_af = '' OR gnomad_af < 0.001)
ORDER BY gene_category, gene_symbol, sample, pos
"""

CATEGORY_ORDER = ["PI3K", "RAS/MAPK", "upstream RTK", "subtype context"]


def read_params(paths):
    out = {}
    for p in paths:
        d = {}
        with open(p) as fh:
            for line in fh:
                k, _, v = line.partition(":")
                d[k.strip()] = v.strip().strip('"')
        sample = re.sub(r"\.somatic_params\.yaml$", "", p.split("/")[-1])
        out[sample] = d
    return out


def protein(hgvsp, hgvsc):
    """ENSP...:p.(Arg88Gln) -> p.Arg88Gln; fall back to the c. change."""
    if hgvsp:
        return hgvsp.split(":")[-1].replace("(", "").replace(")", "")
    return hgvsc.split(":")[-1] if hgvsc else ""


def interpret(r):
    if r["filter"] != "PASS":
        return f"not reported: failed DRAGEN filter ({r['filter']})"
    if r["germline_status"] == "Somatic_Hotspot":
        return "somatic (COSMIC hotspot position)"
    if r["variant_class"] == "somatic":
        return "somatic"
    if r["germline_status"] == "Germline_Proxi":
        return "possible germline (VAF ~0.5, no matched normal)"
    return r["variant_class"]


def short_clinvar(sig):
    """First ClinVar term, abbreviated, for the one-page table."""
    if not sig:
        return "-"
    terms = []
    for t in re.split(r"[;|]", sig):
        t = t.strip().replace("Conflicting classifications of pathogenicity", "Conflicting")
        t = t.replace("Uncertain significance", "VUS")
        if t and t not in terms:
            terms.append(t)
    return "; ".join(terms)


def short_interp(r):
    if r["filter"] != "PASS":
        return f"filtered by DRAGEN ({r['filter'].replace('_', ' ')})"
    if r["germline_status"] == "Somatic_Hotspot":
        return "somatic, hotspot position"
    if r["germline_status"] == "Germline_Proxi":
        return "possibly germline"
    return r["interpretation"]


def fmt_num(x, nd=2):
    try:
        return f"{float(x) + 1e-9:.{nd}f}"
    except (TypeError, ValueError):
        return ""


def main():
    db, out_xlsx, out_md, *param_files = sys.argv[1:]
    params = read_params(param_files)
    con = sqlite3.connect(db)
    con.row_factory = sqlite3.Row

    samples = sorted(params) or [
        r[0] for r in con.execute("SELECT DISTINCT sample FROM panel_variants")
    ]
    genes = list(
        con.execute("SELECT gene_symbol, category, priority, mechanism FROM panel_genes")
    )
    genes.sort(
        key=lambda g: (
            CATEGORY_ORDER.index(g["category"])
            if g["category"] in CATEGORY_ORDER
            else len(CATEGORY_ORDER),
        )
    )
    cand = [dict(r) for r in con.execute(CANDIDATES)]
    for r in cand:
        r["protein"] = protein(r["hgvsp"], r["hgvsc"])
        r["interpretation"] = interpret(r)
    passing = [r for r in cand if r["filter"] == "PASS"]

    # copy-number: any PASS, non-REF segment over the gene
    cnv = {}
    for r in con.execute(
        "SELECT gene_symbol, sample, cnv_call, cnv_filter, copy_ratio FROM panel_gene_cnv"
    ):
        if r["cnv_filter"] == "PASS" and r["cnv_call"] not in ("REF", "."):
            cnv.setdefault((r["gene_symbol"], r["sample"]), set()).add(
                f"copy-number {r['cnv_call']} (ratio {fmt_num(r['copy_ratio'])})"
            )
    n_cnv_pass = len(cnv)

    # ---- gene x sample matrix ----
    cells = {}
    for r in passing:
        txt = f"{r['protein']} (VAF {fmt_num(r['vaf'])})"
        if r["interpretation"].startswith("possible germline"):
            txt += ", possible germline"
        cells.setdefault((r["gene_symbol"], r["sample"]), []).append(txt)
    for k, v in cnv.items():
        cells.setdefault(k, []).extend(sorted(v))

    wb = Workbook()
    bold = Font(bold=True)
    hit_fill = PatternFill("solid", fgColor="FFF2CC")
    wrap = Alignment(wrap_text=True, vertical="top")

    ws = wb.active
    ws.title = "Gene matrix"
    ws.append(["Gene", "Pathway", "Role", "Panel tier"] + samples)
    for c in ws[1]:
        c.font = bold
    for g in genes:
        row = [g["gene_symbol"], g["category"], g["mechanism"], g["priority"]]
        for s in samples:
            row.append("; ".join(cells.get((g["gene_symbol"], s), [])) or "no alteration")
        ws.append(row)
        for c in ws[ws.max_row][4:]:
            if c.value != "no alteration":
                c.fill = hit_fill
    ws.append([])
    ws.append(
        [
            "Protein-altering variants passing all filters (see Methods). "
            + (
                "No copy-number change passed DRAGEN filters in any panel gene."
                if n_cnv_pass == 0
                else "Copy-number changes shown where they passed DRAGEN filters."
            )
        ]
    )
    for col, w in zip("ABCDEFGH", [10, 16, 16, 11] + [34] * len(samples)):
        ws.column_dimensions[col].width = w
    ws.freeze_panes = "B2"

    ws = wb.create_sheet("Alterations")
    head = [
        "Gene", "Sample", "Protein change", "cDNA change", "Consequence",
        "Position (hg38)", "VAF", "Depth", "gnomAD AF", "ClinVar", "Interpretation",
    ]
    ws.append(head)
    for c in ws[1]:
        c.font = bold
    for r in cand:
        ws.append(
            [
                r["gene_symbol"], r["sample"], r["protein"],
                (r["hgvsc"] or "").split(":")[-1], r["consequence"],
                f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}",
                round(float(r["vaf"]), 3) if r["vaf"] not in (None, "") else None,
                r["dp"],
                r["gnomad_af"] if r["gnomad_af"] not in (None, "") else "absent",
                r["clinvar_significance"] or "",
                r["interpretation"],
            ]
        )
    for col, w in zip("ABCDEFGHIJK", [9, 8, 22, 18, 20, 30, 7, 7, 10, 28, 44]):
        ws.column_dimensions[col].width = w

    # ---- methods / limitations text (shared by xlsx and md) ----
    cov = ", ".join(f"{s} {float(params[s]['median_cov']):.0f}x" for s in samples if s in params)
    thr = ", ".join(f"{s} {params[s]['callability_thresh']}" for s in samples if s in params)
    width = {params[s].get("cnv_interval_width") for s in samples if s in params}
    width = "/".join(sorted(w for w in width if w)) or "1000"
    width_kb = f"{int(width) // 1000} kb" if width.isdigit() else f"{width} bp"
    n_genes = len(genes)

    methods = (
        "Whole-genome sequencing data from organoids "
        f"{', '.join(samples[:-1]) + ' and ' + samples[-1] if len(samples) > 1 else samples[0]} "
        f"(median autosomal coverage {cov}) were analysed with Illumina DRAGEN 4.4.7 "
        "in somatic tumor-only mode, as no matched normal tissue was available. "
        "Reads were adapter- and quality-trimmed, aligned to the hg38 reference "
        "(alt-masked, linear) and duplicate-marked. Small variants were called with "
        "the DRAGEN somatic caller using Illumina's WGS systematic-noise baseline; "
        "variants present in gnomAD or 1000 Genomes were tagged as germline. "
        f"The minimum depth for a callable position was set to half the median "
        f"coverage ({thr} reads), and COSMIC hotspot positions were exempt from "
        f"germline tagging. Copy number was called with self-normalization in "
        f"{width_kb} bins. Variants were annotated with Illumina Connected "
        "Annotations 3.26 (Ensembl canonical transcripts, gnomAD, ClinVar). We "
        f"report variants in {n_genes} genes of the PI3K and RAS/MAPK pathways, "
        "upstream receptor tyrosine kinases, and endometrial cancer molecular-subtype "
        "genes that passed all DRAGEN filters, altered the protein sequence (high or "
        "moderate predicted impact), had a gnomAD allele frequency below 0.1% and "
        "were not tagged germline by population databases."
    )
    limits = [
        "No matched normal: somatic status is inferred, not proven. A rare variant "
        "near 50% VAF cannot be told apart from a heterozygous germline variant in "
        "a near-pure organoid.",
        "Coverage ~25x: variants below ~15% VAF have 4 or fewer supporting reads "
        "and may be missed.",
        "Copy number: no event passed DRAGEN filters and tumor purity could not be "
        f"estimated; focal changes spanning only a few {width_kb} bins may be missed.",
    ]

    ws = wb.create_sheet("Methods")
    ws.column_dimensions["A"].width = 110
    for line in ["Methods", methods, "", "Limitations"] + [f"- {x}" for x in limits]:
        ws.append([line])
        ws.cell(ws.max_row, 1).alignment = wrap
    ws["A1"].font = bold
    ws["A4"].font = bold
    wb.save(out_xlsx)

    # ---- one-page Markdown ----
    by_sample = {s: [r for r in passing if r["sample"] == s] for s in samples}
    hit_genes = {r["gene_symbol"] for r in passing}
    core_clean = [
        g["gene_symbol"]
        for g in genes
        if g["priority"] == "core" and g["gene_symbol"] not in hit_genes
    ]
    lines = [
        "# PI3K / RAS-MAPK alterations in endometrial cancer organoids",
        "",
        "## Result",
        "",
    ]
    for s in samples:
        som = [
            f"{r['gene_symbol']} {r['protein']}"
            for r in by_sample[s]
            if r["interpretation"].startswith("somatic")
        ]
        unres = [
            f"{r['gene_symbol']} {r['protein']}"
            for r in by_sample[s]
            if r["interpretation"].startswith("possible germline")
        ]
        txt = f"- **{s}**: " + (", ".join(som) if som else "no somatic alteration")
        if unres:
            txt += f"; {', '.join(unres)} (possibly germline)"
        lines.append(txt + ".")
    lines += [
        f"- No alteration in the other core genes ({', '.join(core_clean)}).",
        "- No copy-number change in any panel gene."
        if n_cnv_pass == 0
        else "- Copy-number changes: see the gene matrix.",
        "",
        "| Gene | Sample | Protein change | VAF | Depth | ClinVar | Interpretation |",
        "|-------|--------|-------------------|-----|------|----------------|-----------------------|",
    ]
    for r in cand:
        lines.append(
            f"| {r['gene_symbol']} | {r['sample']} | {r['protein']} | "
            f"{fmt_num(r['vaf'])} | {r['dp']} | "
            f"{short_clinvar(r['clinvar_significance'])} | "
            f"{short_interp(r)} |"
        )
    lines += ["", "## Methods", "", methods, "", "## Limitations", ""]
    lines += [f"- {x}" for x in limits]
    lines += [
        "",
        f"Full per-gene results for all {n_genes} genes: "
        f"`{out_xlsx.split('/')[-1]}`.",
        "",
    ]
    with open(out_md, "w") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    main()
