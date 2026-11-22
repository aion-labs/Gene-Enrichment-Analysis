"""Supplementary Table 8 (HTML): sensitivity of iGEA results to the analysis settings.
Input: 07 Reviewer Analyses/step10_extras/s10a_parameter_sensitivity.tsv. Usage: python build_supp_table_15.py <tsv> <out.html>"""
import sys, html
import pandas as pd

d = pd.read_csv(sys.argv[1], sep="\t")
DESC = {"default": ("Default", "p < 0.01 · overlap ≥ 3 · 10–600 genes · ≤ 30 iterations"),
        "p < 0.005": ("Stricter selection", "p < 0.005"), "p < 0.05": ("Looser selection", "p < 0.05"),
        "min overlap 2": ("Smaller minimum overlap", "≥ 2 input genes per term"),
        "min overlap 5": ("Larger minimum overlap", "≥ 5 input genes per term"),
        "set size 15-500": ("Narrower gene-set sizes", "15–500 genes"), "set size 5-1000": ("Wider gene-set sizes", "5–1,000 genes"),
        "max 10 iterations": ("Fewer iterations", "≤ 10 per collection")}
GROUP = {"default": "", "p < 0.005": "Selection threshold", "p < 0.05": "Selection threshold", "min overlap 2": "Minimum overlap",
         "min overlap 5": "Minimum overlap", "set size 15-500": "Gene-set size range", "set size 5-1000": "Gene-set size range",
         "max 10 iterations": "Iteration cap"}

def pbadge(p, kept):
    if kept == 0:
        return '<span class="pb na">no network</span>'
    if p < 0.001:
        return '<span class="pb sig">&lt; 0.001</span>'
    return f'<span class="pb {"sig" if p < 0.05 else "ns"}">{p:.2f}</span>' if p >= 0.1 else f'<span class="pb {"sig" if p < 0.05 else "ns"}">{p:.3f}</span>'

def bar(v, kept):
    if kept == 0:
        return '<span class="muted">–</span>'
    w = max(2, round(100 * v))
    return f'<div class="bar"><i style="width:{w}%"></i><span>{v:.2f}</span></div>'

def table(ds, sig_col_note):
    rows = []
    for r in d[d.dataset == ds].itertuples():
        name, val = DESC[r.setting]
        cls = ' class="def"' if r.setting == "default" else ""
        rows.append(f"<tr{cls}><td><b>{html.escape(name)}</b><div class='sub'>{html.escape(val)}</div></td>"
                    f"<td class='n'>{r.selected}</td><td class='n'><b>{r.kept}</b></td><td class='n'>{r.chance_terms_per_list:.1f}</td>"
                    f"<td class='n'>{r.cluster_genes if r.kept else '–'}</td><td class='n'>{r.cluster_terms if r.kept else '–'}</td>"
                    f"<td>{pbadge(r.p_cluster_size, r.kept)}</td><td>{pbadge(r.p_edges_per_gene, r.kept)}</td>"
                    f"<td>{bar(r.jaccard_kept_terms_vs_default, r.kept)}</td><td>{bar(r.jaccard_genes_vs_default, r.kept)}</td></tr>")
    return ("<div class='tbl'><table><thead><tr><th>Setting</th><th class='n'>Terms selected</th><th class='n'>Terms kept<br><span>q &lt; 0.10</span></th>"
            "<th class='n'>Chance terms<br><span>per random list</span></th><th class='n'>Cluster<br><span>genes</span></th>"
            "<th class='n'>Cluster<br><span>terms</span></th><th>Network p<br><span>cluster size</span></th><th>Network p<br><span>edges per gene</span></th>"
            "<th>Same terms as default<br><span>Jaccard, kept terms</span></th><th>Same genes as default<br><span>Jaccard, genes</span></th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table></div>")

hiv = d[d.dataset == "HIV"]; rnd = d[d.dataset == "Random control"]
rng = lambda s: f"{s.min():.2f}–{s.max():.2f}"
nd = hiv[hiv.setting != "default"]
CSS = """:root{--ink:#22252b;--muted:#5d6470;--line:#e3e6ea;--bg:#fbfbfa;--a:#2a78d6;--a-bg:#eef5fd;--b:#1baf7a;--b-bg:#ecf8f3;--c:#eb6834;--c-bg:#fdf1ec}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
main{max-width:1180px;margin:0 auto;padding:40px 28px 60px}.kicker{font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:0}
h1{font-size:26px;margin:4px 0 12px}h2{font-size:19px;margin:34px 0 4px;padding-bottom:6px;border-bottom:2px solid var(--line)}
p{margin:8px 0}.lead{color:#3a3f47;max-width:900px}.muted{color:var(--muted)}
.cards{display:flex;gap:12px;flex-wrap:wrap;margin:22px 0 8px}.card{flex:1 1 260px;background:#fff;border:1px solid var(--line);border-radius:10px;padding:12px 15px}
.card b{display:block;font-size:14px;margin-bottom:3px}.card span{font-size:13px;color:var(--muted)}
.card.a{border-top:4px solid var(--a)}.card.b{border-top:4px solid var(--b)}.card.c{border-top:4px solid var(--c)}
.tbl{overflow-x:auto;margin:12px 0}table{border-collapse:collapse;width:100%;background:#fff;border:1px solid var(--line);font-size:13.5px}
th{background:#f3f4f6;text-align:left;font-weight:600;padding:8px 10px;border-bottom:1px solid var(--line);vertical-align:bottom;font-size:12.5px}
th span{font-weight:400;color:var(--muted);font-size:11.5px}th.n,td.n{text-align:right}
td{padding:8px 10px;border-bottom:1px solid var(--line);vertical-align:middle}tr:last-child td{border-bottom:none}
tr.def td{background:#f6f9fe}tr.def td:first-child{box-shadow:inset 4px 0 0 var(--a)}.sub{font-size:12px;color:var(--muted)}
.pb{display:inline-block;font-size:12px;font-weight:600;border-radius:999px;padding:1px 9px;white-space:nowrap}
.pb.sig{background:var(--b-bg);color:#10845a;border:1px solid #bfe6d6}.pb.ns{background:#f3f4f6;color:#5d6470;border:1px solid var(--line)}
.pb.na{background:#fff;color:#9aa1ab;border:1px dashed #d5d9de;font-weight:500}
.bar{position:relative;width:130px;height:18px;background:#f1f3f5;border-radius:4px;overflow:hidden}.bar i{position:absolute;left:0;top:0;bottom:0;background:#9cc3ef}
.bar span{position:relative;font-size:12px;font-weight:600;padding-left:6px;line-height:18px}
.notes{font-size:13px;color:#3a3f47;max-width:960px}.notes li{margin:4px 0}
footer{margin-top:30px;padding-top:12px;border-top:1px solid var(--line);color:var(--muted);font-size:12.5px}
@media print{body{background:#fff;font-size:10.5pt}main{padding:0;max-width:none}table{break-inside:avoid}h2{break-after:avoid}}"""
doc = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Supplementary Table 8 — Sensitivity to analysis settings</title><style>{CSS}</style></head><body><main>
<p class="kicker">Supplementary Table 8</p><h1>Sensitivity of iGEA results to the analysis settings</h1>
<p class="lead">Each row changes one setting from the default analysis. For every setting, the random gene lists were generated and analysed
with the same changed setting, so the term-level q-values and the network test always match the analysis that was run.</p>
<div class="cards">
<div class="card b"><b>The HIV network stays significant</b><span>Empirical p &lt; 0.001 for the largest cluster under every setting
(edges per gene p {nd.p_edges_per_gene.min():.3f}–{nd.p_edges_per_gene.max():.3f}).</span></div>
<div class="card c"><b>The random list never becomes significant</b><span>Network p ≥ {rnd[rnd.kept > 0].p_cluster_size.min():.2f}; at most
{rnd.kept.max()} term kept under any setting.</span></div>
<div class="card a"><b>Settings change terms more than genes</b><span>Overlap with the default: kept terms {rng(nd.jaccard_kept_terms_vs_default)},
genes {rng(nd.jaccard_genes_vs_default)}.</span></div></div>
<h2>HIV dependency factors (265 genes)</h2>{table("HIV", "")}
<h2>Random control (264 genes drawn at random from the background)</h2>
<p class="muted" style="font-size:13.5px">The negative control of the Results, analysed in exactly the same way. "No network": no term passed the term-level filter.</p>
{table("Random control", "")}
<h2>Notes</h2><ul class="notes">
<li><b>Default analysis:</b> one-sided Fisher's exact test; p &lt; 0.01 for selection; at least 3 input genes per term; gene sets of 10–600 genes
after restriction to the background (23,391 genes); at most 30 iterations per collection; the ten collections of Supplementary Table 2.</li>
<li><b>Random gene lists:</b> 2,000 lists of the input size and 200 per reduced size for the q-value adjustment, generated with the setting of each
row (seed 2026). With 2,000 lists the smallest attainable empirical p = (r + 1)/(n + 1) is 0.0005, shown as &lt; 0.001. The analyses in the main text use 10,000 lists.</li>
<li><b>Terms kept:</b> terms with empirical q &lt; 0.10, which form the network. <b>Chance terms per random list:</b> mean number of terms selected
from a random list of the same size with the same setting.</li>
<li><b>Same terms / genes as default:</b> Jaccard similarity between the kept terms (or the genes assigned to them) of each setting and of the default analysis of the same gene list.</li>
</ul>
<footer>iGEA, per-run statistics (per_run_null.run_statistics). Script: 07 Reviewer Analyses/step10_extras/s10a_sensitivity.py; data: s10a_parameter_sensitivity.tsv.</footer>
</main></body></html>"""
open(sys.argv[2], "w").write(doc)
print("written", sys.argv[2])
