"""Supplementary Table 2: the ten gene-set collections (HTML + tab-separated text).
Input: collection_table_data.tsv (from collection_table_data.py: current code, one universe, one-sided Fisher;
HIV counts read from Supplementary Tables 4 and 5)."""
import csv, html, os, re, sys
SRC, OUT = sys.argv[1], sys.argv[2]   # collection_table_data.tsv, output folder
rows = {r["Collection"]: r for r in csv.DictReader(open(SRC), delimiter="\t")}
ORDER = [("H: Hallmark Gene Sets", "Hallmark", "H"),
         ("C2: BioCarta", "BioCarta", "C2:CP"), ("C2: KEGG Legacy", "KEGG Legacy", "C2:CP"),
         ("C2: KEGG MEDICUS", "KEGG MEDICUS", "C2:CP"), ("C2: Pathway Interaction Database", "Pathway Interaction Database", "C2:CP"),
         ("C2: Reactome Pathways", "Reactome", "C2:CP"), ("C2: WikiPathways", "WikiPathways", "C2:CP"),
         ("C5: Gene Ontology: Biological Process", "GO Biological Process", "C5:GO"),
         ("C5: Gene Ontology: Cellular Component", "GO Cellular Component", "C5:GO"),
         ("C5: Gene Ontology: Molecular Function", "GO Molecular Function", "C5:GO")]
COLS = [("Gene sets in file", "Gene sets in file"), ("Gene sets used (10-600 genes)", "Gene sets used (10–600 genes)"),
        ("Genes covered", "Genes covered"), ("Median set size", "Median set size")]
HIV = [("HIV genes covered", "Input genes covered"), ("HIV: possible terms (>= 3 genes)", "Possible terms (≥ 3 genes)"),
       ("HIV: ORA p < 0.01", "ORA p < 0.01"), ("HIV: ORA FDR < 0.05", "ORA FDR < 0.05"),
       ("HIV: iGEA selected", "iGEA selected"), ("HIV: iGEA kept (q < 0.10)", "iGEA kept (q < 0.10)")]
UNION_GENES, HIV_ANY, HIV_N = "18,506", "246", "265"
fmt = lambda v: f"{int(v):,}"
tot = {k: sum(int(rows[c][k]) for c, _, _ in ORDER) for k, _ in COLS + HIV}
total_cells = {"Gene sets in file": fmt(tot["Gene sets in file"]), "Gene sets used (10-600 genes)": fmt(tot["Gene sets used (10-600 genes)"]),
               "Genes covered": UNION_GENES + "*", "Median set size": "—", "HIV genes covered": f"{HIV_ANY}* of {HIV_N}"}
for k, _ in HIV[1:]:
    total_cells[k] = fmt(tot[k])

# ---------------------------------------------------------------- text version
with open(f"{OUT}/Supp Table 02 - Gene-set collections.txt", "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["Collection", "MSigDB category", "GMT file (MSigDB v2025.1)"] + [h for _, h in COLS] + ["HIV: " + h for _, h in HIV])
    for c, name, cat in ORDER:
        r = rows[c]
        w.writerow([name, cat, r["File"]] + [r[k] for k, _ in COLS + HIV])
    w.writerow(["Total (10 collections)", "", ""] + [total_cells[k].replace(",", "") for k, _ in COLS + HIV])
    w.writerow([])
    for note in ["* Number of distinct genes, not the sum over collections.",
                 "Gene sets used: restricted to the background (23,391 validated symbols), then 10-600 genes. Genes covered = collection-specific background size M.",
                 "Possible terms: used gene sets with at least 3 HIV genes. ORA: one-sided Fisher's exact test, p < 0.01, overlap >= 3, BH FDR per collection (Supplementary Table 4).",
                 "iGEA: selected at p < 0.01; kept at empirical q < 0.10 (10,000 random gene lists); kept terms form the network (Supplementary Table 5)."]:
        w.writerow([note])

# ---------------------------------------------------------------- HTML
css = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "st02.css")).read()
css += """
.tbl{overflow-x:auto;margin:16px 0}
main{max-width:1240px} table.ct{font-size:13px} table.ct td,table.ct th{padding:6px 8px} table.ct td:first-child{white-space:nowrap}
table.ct th{font-size:12px;vertical-align:bottom;line-height:1.25;max-width:92px}
table.ct th.grp{text-align:center;font-size:12px;letter-spacing:.04em;text-transform:uppercase;color:#5d6470;background:#fff;border-bottom:none}
table.ct th.grp.hiv{background:#eef5fd;color:#2a78d6}
table.ct td.n,table.ct th.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap} table.ct th.n{white-space:normal}
table.ct td.hiv,table.ct th.hiv{background:#f6f9fe}
table.ct tr.sec td{background:#f3f4f6;font-size:11.5px;font-weight:600;letter-spacing:.05em;text-transform:uppercase;color:#5d6470;padding:5px 12px}
table.ct tr.tot td{font-weight:700;border-top:2px solid #c9ced6;background:#fafafa}
table.ct tr.tot td.hiv{background:#eef5fd}
.cat{display:inline-block;font-size:11px;font-weight:600;border-radius:4px;padding:1px 6px;margin-left:6px;color:#5d6470;background:#eef0f3}
.notes{font-size:13.5px;color:#3a3f47;margin-top:18px}
.notes li{margin:5px 0}
"""
SECTIONS = {"H": "Hallmark", "C2:CP": "C2 · Canonical pathways (six subcollections)", "C5:GO": "C5 · Gene Ontology"}
body, last = [], None
for c, name, cat in ORDER:
    if cat != last:
        body.append(f'<tr class="sec"><td colspan="11">{SECTIONS[cat]}</td></tr>')
        last = cat
    r = rows[c]
    cells = "".join(f'<td class="n">{fmt(r[k]) if k != "Median set size" else r[k]}</td>' for k, _ in COLS)
    cells += "".join(f'<td class="n hiv">{fmt(r[k])}</td>' for k, _ in HIV)
    body.append(f"<tr><td>{html.escape(name)}</td>{cells}</tr>")
tcells = "".join(f'<td class="n">{total_cells[k]}</td>' for k, _ in COLS) + \
         "".join(f'<td class="n hiv">{total_cells[k]}</td>' for k, _ in HIV)
body.append(f'<tr class="tot"><td>Total (10 collections)</td>{tcells}</tr>')
head1 = '<tr><th class="grp" colspan="5">Collection (MSigDB v2025.1)</th><th class="grp hiv" colspan="6">HIV dependency factors (265 genes)</th></tr>'
head2 = "<tr><th>Collection</th>" + "".join(f'<th class="n">{h}</th>' for _, h in COLS) + \
        "".join(f'<th class="n hiv">{h}</th>' for _, h in HIV) + "</tr>"
page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Supplementary Table 2 — Gene-set collections</title><style>{css}</style></head>
<body><main>
<p class="kicker">Supplementary Table 2</p>
<h1>Gene-set collections</h1>
<p class="lead">The ten MSigDB collections used in every analysis, figure and set of random gene lists in this study, unless a subset is stated. The aggregate Canonical Pathways collection (C2:CP), which contains the six canonical pathway subcollections, was not used, so that no gene set enters the analysis twice. The right-hand columns show, for the HIV dependency factors, how many terms each collection could select and how many ORA and iGEA select.</p>
<div class="tbl"><table class="ct"><thead>{head1}{head2}</thead><tbody>{''.join(body)}</tbody></table></div>
<ul class="notes">
<li><b>Gene sets used:</b> each gene set is first restricted to the background (23,391 validated human gene symbols) and kept if it then has 10–600 genes.</li>
<li><b>Genes covered:</b> background genes in at least one used gene set; this is the collection-specific background size M of the enrichment test. *Totals of genes covered are the number of distinct genes, not the sum.</li>
<li><b>Input genes covered:</b> HIV genes (265 in the background) present in the collection; this is the collection-specific input size N at the first iteration.</li>
<li><b>Possible terms:</b> used gene sets that contain at least three HIV genes (the minimum overlap), i.e. the largest number of terms that could reach significance.</li>
<li><b>ORA:</b> one-sided Fisher's exact test, p &lt; 0.01, minimum overlap 3; FDR = Benjamini–Hochberg within each collection (Supplementary Table 4).</li>
<li><b>iGEA:</b> terms selected by the iterative procedure (p &lt; 0.01) and kept by the term-level filter (empirical q &lt; 0.10, 10,000 random gene lists); kept terms form the network (Supplementary Table 5).</li>
</ul>
<footer>Gene sets: MSigDB v2025.1, human symbols harmonized to current NCBI gene symbols. Counts computed with the iGEA code (iGEA repository, publication/csbj_2026/build/supplement/collection_table_data.py and build_supp_table_02.py).</footer>
</main></body></html>"""
open(f"{OUT}/Supp Table 02 - Gene-set collections.html", "w", encoding="utf8").write(page)
print("ok")
