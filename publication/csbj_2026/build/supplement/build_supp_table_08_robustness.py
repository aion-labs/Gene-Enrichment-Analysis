"""Supplementary Table 8 (HTML + TXT): robustness of iGEA results.
Part A: analysis settings (step10_extras/s10a_parameter_sensitivity.tsv).
Part B: input subsampling and near ties (step06_stability/s06_subsamples.tsv, s06_term_recurrence.tsv, s06_near_ties.tsv).
Part C: redundancy between collections (step03_redundancy/s03_redundancy_matrix_J90.csv, s03_dedup_summary.json,
        s03_sensitivity_dedup.csv, s03_sensitivity_loo.csv).
Part D (R2, Dec 2026): annotation-matched random gene lists (step10_extras/s10c_matched_*.tsv, s10c_matched_HIV_terms.tsv),
        compared with the standard random lists of the main analyses (step09_end_to_end).
Usage: python build_supp_table_08_robustness.py <analyses dir> <out dir>"""
import sys, os, json, html
import pandas as pd
A, OUT = sys.argv[1], sys.argv[2]
CSS = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "st08.css")).read()
d = pd.read_csv(f"{A}/step10_extras/s10a_parameter_sensitivity.tsv", sep="\t")
sub = pd.read_csv(f"{A}/step06_stability/s06_subsamples.tsv", sep="\t")
rec = pd.read_csv(f"{A}/step06_stability/s06_term_recurrence.tsv", sep="\t")
nt = pd.read_csv(f"{A}/step06_stability/s06_near_ties.tsv", sep="\t")
mat = pd.read_csv(f"{A}/step03_redundancy/s03_redundancy_matrix_J90.csv")
dd = json.load(open(f"{A}/step03_redundancy/s03_dedup_summary.json"))
dedup = pd.read_csv(f"{A}/step03_redundancy/s03_sensitivity_dedup.csv")
loo = pd.read_csv(f"{A}/step03_redundancy/s03_sensitivity_loo.csv")
E = html.escape
import glob, re
# ---------------------------------------------------------------- Part D: annotation-matched random lists
RUNS = {"HIV": "HIV_dependency_factors.symbols", "Hypoxia": "hypoxia-genes.symbols", "Vaccination": "vaccination_NAKAYA_HAI_28DY_3DY_POSITIVE"}
drows, dtxt = [], []
for ds, pre in RUNS.items():
    mt = pd.read_csv(f"{A}/step10_extras/s10c_matched_{ds}.tsv", sep="\t").set_index("Metric")
    run = glob.glob(f"{A}/step09_end_to_end/{pre}_2026*/")[0]
    st = pd.read_csv(glob.glob(f"{run}*_statistics.tsv")[0], sep="\t").set_index("Metric")
    rep = open(glob.glob(f"{run}*_statistical_report.txt")[0]).read()
    kept_u, sel = map(int, re.search(r"Terms kept \(q < 0.1\): (\d+) of (\d+)", rep).groups())
    ch_u = float(re.search(r"selects on average ([\d.]+) terms", rep).group(1))
    m0 = mt.iloc[0]
    r = dict(ds=ds, sel=sel, ann_in=m0.input_mean_annotation, ann_bg=m0.background_mean_annotation, ch_u=ch_u, ch_m=m0.chance_terms_per_list,
             kept_u=kept_u, kept_m=int(m0.terms_kept), p_cs_u=st.loc["Cluster Size", "Empirical p"], p_cs_m=mt.loc["Cluster Size", "Empirical p"],
             p_epg_u=st.loc["Average Edges per Gene", "Empirical p"], p_epg_m=mt.loc["Average Edges per Gene", "Empirical p"],
             genes_m=int(mt.loc["Number of Genes", "Observed"]), terms_m=int(mt.loc["Number of Terms", "Observed"]))
    drows.append(r)
def pd_(p):
    if p < 0.0006: return '<span class="pb sig">&lt; 0.001</span>'
    return f'<span class="pb {"sig" if p < 0.05 else "ns"}">{p:.3f}</span>' 
DROWS = "".join(f"<tr><td><b>{r['ds']}</b></td><td class='n'>{r['ann_in']:.1f}</td><td class='n'>{r['ch_u']:.1f} → {r['ch_m']:.1f}</td>"
                f"<td class='n'>{r['kept_u']} → <b>{r['kept_m']}</b> <span class='muted'>of {r['sel']}</span></td><td class='n'>{r['genes_m']} / {r['terms_m']}</td>"
                f"<td>{pd_(r['p_cs_u'])} → {pd_(r['p_cs_m'])}</td><td>{pd_(r['p_epg_u'])} → {pd_(r['p_epg_m'])}</td></tr>" for r in drows)
mterms = pd.read_csv(f"{A}/step10_extras/s10c_matched_HIV_terms.tsv", sep="\t")
hrun = glob.glob(f"{A}/step09_end_to_end/HIV_dependency_factors.symbols_2026*/")[0]
hq = pd.read_csv(glob.glob(f"{hrun}*_term_qvalues.tsv")[0], sep="\t")
hq["Community"] = pd.read_csv(f"{hrun}combined_iterative_results.tsv", sep="\t")["Community"].values
hx = hq.merge(mterms, on=["Library", "Term", "Iteration"], suffixes=("", "_m")); hx = hx[hx.Kept]
assert len(hx) == 56
un = pd.read_csv(f"{A}/step13_r2/s13a_terms_HIV.tsv", sep="\t").set_index(["Library", "Term"])["q-value (no signal adjustment)"]
hx["q_unadj"] = [un[(a, b)] for a, b in zip(hx.Library, hx.Term)]
dropped = hx[~hx.Kept_m]; n_unadj = int((hx.q_unadj >= 0.10).sum()); n_both = int(((hx.q_unadj >= 0.10) & ~hx.Kept_m).sum())
CNAME = {"C1": "Golgi vesicle transport", "C2": "growth-factor and NF-κB signaling, endocytosis", "C3": "Mediator-dependent transcription",
         "C4": "nuclear pore, nuclear export", "C5": "RNA processing, RNA exosome", "C6": "ER glycosylation, GPI anchor", "C7": "cullin-RING ubiquitin ligases"}
cg = hx.groupby("Community").agg(n=("Term", "size"), m=("Kept_m", "sum"))
CROWS = "".join(f"<tr><td><b>{c}</b> {E(CNAME[c])}</td><td class='n'>{int(r.n)}</td><td class='n'><b>{int(r.m)}</b></td></tr>" for c, r in cg.iterrows())
# ---------------------------------------------------------------- Part A
DESC = {"default": ("Default", "p < 0.01; overlap ≥ 3; 10–600 genes; ≤ 30 iterations"),
        "p < 0.005": ("Stricter selection", "p < 0.005"), "p < 0.05": ("Looser selection", "p < 0.05"),
        "min overlap 2": ("Smaller minimum overlap", "≥ 2 input genes per term"), "min overlap 5": ("Larger minimum overlap", "≥ 5 input genes per term"),
        "set size 15-500": ("Narrower gene-set sizes", "15–500 genes"), "set size 5-1000": ("Wider gene-set sizes", "5–1,000 genes"),
        "max 10 iterations": ("Fewer iterations", "≤ 10 per collection")}
def pb(p, ok=True):
    if not ok: return '<span class="pb na">no network</span>'
    if p < 0.0006: return '<span class="pb sig">&lt; 0.001</span>'
    return f'<span class="pb {"sig" if p < 0.05 else "ns"}">{p:.3f}</span>' if p < 0.1 else f'<span class="pb ns">{p:.2f}</span>'
def bar(v, ok=True):
    if not ok: return '<span class="muted">–</span>'
    return f'<div class="bar"><i style="width:{max(2, round(100 * v))}%"></i><span>{v:.2f}</span></div>'
def tableA(ds):
    rows = []
    for r in d[d.dataset == ds].itertuples():
        nm, val = DESC[r.setting]; ok = r.kept > 0
        rows.append(f"<tr{' class=def' if r.setting == 'default' else ''}><td><b>{E(nm)}</b><div class='sub'>{E(val)}</div></td>"
                    f"<td class='n'>{r.selected}</td><td class='n'><b>{r.kept}</b></td><td class='n'>{r.chance_terms_per_list:.1f}</td>"
                    f"<td class='n'>{r.cluster_genes if ok else '–'}</td><td class='n'>{r.cluster_terms if ok else '–'}</td>"
                    f"<td>{pb(r.p_cluster_size, ok)}</td><td>{pb(r.p_edges_per_gene, ok)}</td>"
                    f"<td>{bar(r.jaccard_kept_terms_vs_default, ok)}</td><td>{bar(r.jaccard_genes_vs_default, ok)}</td></tr>")
    return ("<div class='tbl'><table><thead><tr><th>Setting</th><th class='n'>Terms selected</th><th class='n'>Terms kept<br><span>q &lt; 0.10</span></th>"
            "<th class='n'>Chance terms<br><span>per random list</span></th><th class='n'>Cluster<br><span>genes</span></th><th class='n'>Cluster<br><span>terms</span></th>"
            "<th>Network p<br><span>cluster size</span></th><th>Network p<br><span>edges per gene</span></th><th>Same terms as default<br><span>Jaccard, kept terms</span></th>"
            f"<th>Same genes as default<br><span>Jaccard, genes</span></th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>")
hiv = d[d.dataset == "HIV"]; rnd = d[d.dataset == "Random control"]; nd = hiv[hiv.setting != "default"]
# ---------------------------------------------------------------- Part B
def mr(c, f="{:.0f}"): return f"{f.format(sub[c].median())} ({f.format(sub[c].min())}–{f.format(sub[c].max())})"
k = rec[rec.kept_full]; ntt = nt["near tie (<=2x)"].astype(bool)
B_ROWS = [("Terms selected (p < 0.01)", mr("selected"), "62"), ("Terms kept (q < 0.10)", mr("kept"), "56"),
          ("Full-list kept terms recovered", mr("full_kept_recovered"), "of 56"), ("Kept terms not kept in the full list", mr("new_kept_terms"), "–"),
          ("Largest cluster: genes", mr("cluster_genes"), "168"), ("Largest cluster: terms", mr("cluster_terms"), "56"),
          ("Cluster genes shared with the full-list cluster (Jaccard)", mr("cluster_gene_jaccard_vs_full", "{:.2f}"), "1"),
          ("Network p, number of genes", f"{sub.network_p_genes.max():.4f} in all {len(sub)}" if sub.network_p_genes.nunique() == 1 else f"≤ {sub.network_p_genes.max():.4f}", "0.0001")]
brows = "".join(f"<tr><td>{E(a)}</td><td class='n'>{E(b)}</td><td class='n'>{E(c)}</td></tr>" for a, b, c in B_ROWS)
low = k[k.recurrence_kept < 0.5].sort_values("recurrence_kept")
lowrows = "".join(f"<tr><td>{E(r.Library.split(': ')[-1])}</td><td>{E(r.Term)}</td><td class='n'>{r.q_full:.3f}</td><td class='n'>{r.recurrence_kept:.2f}</td></tr>" for r in low.itertuples())
# ---------------------------------------------------------------- Part C
cols = [c for c in mat.columns if c not in ("Collection", "Gene sets (10-600)", "Any other (J>=0.9)", "Any other (identical)")]
mhead = "".join(f"<th class='n'>{E(c)}</th>" for c in cols)
mrows = "".join("<tr><td><b>" + E(r["Collection"]) + f"</b></td><td class='n'>{int(r['Gene sets (10-600)']):,}</td>" +
                "".join(f"<td class='n'>{'' if pd.isna(r[c]) or r[c] == '' else int(r[c])}</td>" for c in cols) +
                f"<td class='n'><b>{int(r['Any other (J>=0.9)'])}</b></td><td class='n'>{int(r['Any other (identical)'])}</td></tr>" for _, r in mat.iterrows())
n_any = int(mat["Any other (J>=0.9)"].sum()); n_sets = int(mat["Gene sets (10-600)"].sum()); n_drop = sum(dd["removed_by_dedup"].values())
drop_txt = ", ".join(f"{k_} {v}" for k_, v in dd["removed_by_dedup"].items() if v)
def crow(lab, r):
    return (f"<tr><td>{E(lab)}</td><td class='n'>{int(r['iGEA terms'])}</td><td class='n'><b>{int(r['Kept (q < 0.10)'])}</b></td>"
            f"<td class='n'>{int(r['Cluster genes'])}</td><td class='n'>{int(r['Cluster terms'])}</td><td class='n'>{r['Edges/gene']:.2f}</td>"
            f"<td>{pb(r['p Cluster Size'])}</td><td>{pb(r['p Average Edges per Gene'])}</td></tr>")
crows = crow("All ten, deduplicated", dedup.iloc[0]) + "".join(crow(r["Configuration"][0].upper() + r["Configuration"][1:], r) for _, r in loo.iterrows())
CHEAD = ("<tr><th>Collections</th><th class='n'>Terms selected</th><th class='n'>Terms kept</th><th class='n'>Cluster genes</th><th class='n'>Cluster terms</th>"
         "<th class='n'>Edges per gene</th><th>Network p<br><span>cluster size</span></th><th>Network p<br><span>edges per gene</span></th></tr>")
doc = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Supplementary Table 8 — Robustness of iGEA results</title><style>{CSS}</style></head><body><main>
<p class="kicker">Supplementary Table 8</p><h1>Robustness of iGEA results</h1>
<p class="lead">Four analyses test whether the HIV results depend on choices that a user, or the data, could have made differently:
(A) the analysis settings, (B) the composition of the input list and near ties in the iterative selection, (C) redundancy between
gene-set collections, and (D) the annotation density of the input genes. In every analysis the random gene lists were generated and analyzed in the same way as the gene list being tested.</p>
<div class="cards">
<div class="card b"><b>A · Settings</b><span>The HIV network stays significant under every setting (p &lt; 0.001); the random control never does
(p ≥ {rnd[rnd.kept > 0].p_cluster_size.min():.2f}).</span></div>
<div class="card a"><b>B · Subsampling</b><span>All {len(sub)} subsamples of 90% of the list give a significant network; cluster genes overlap the
full-list cluster with median Jaccard {sub.cluster_gene_jaccard_vs_full.median():.2f}.</span></div>
<div class="card c"><b>C · Collection redundancy</b><span>{n_any} of {n_sets:,} gene sets have a near-duplicate in another collection; removing
duplicates or any one collection leaves the network significant.</span></div>
<div class="card d"><b>D · Annotation density</b><span>With random lists matched to the annotation density of the input, all three networks stay significant,
but fewer terms are kept (HIV {drows[0]['kept_u']} → {drows[0]['kept_m']}).</span></div></div>

<h2>A. Analysis settings</h2>
<p class="muted" style="font-size:13.5px">Each row changes one setting from the default analysis; the random lists use the same changed setting.</p>
<h3>HIV dependency factors (265 genes)</h3>{tableA("HIV")}
<h3>Random control (264 genes drawn at random from the background)</h3>{tableA("Random control")}

<h2>B. Input subsampling and near ties</h2>
<p class="muted" style="font-size:13.5px">Selecting a term and removing its genes changes all later iterations, so the result is path dependent.
{len(sub)} random subsamples of 90% of the HIV list ({int(round(0.9 * 265))} genes; seed 2026) were analyzed with the default settings and compared
with the full-list analysis. All subsamples share one set of 5,000 random lists of their size (the smallest attainable empirical p is 0.0002).</p>
<div class='tbl'><table><thead><tr><th>Quantity</th><th class='n'>Subsamples<br><span>median (range)</span></th><th class='n'>Full list</th></tr></thead><tbody>{brows}</tbody></table></div>
<p class="notes">Of the 56 terms kept in the full list, {int((k.recurrence_kept >= 0.8).sum())} were kept in at least 80% of the subsamples and
{int((k.recurrence_kept >= 0.5).sum())} in at least half (median recurrence {k.recurrence_kept.median():.2f}).
Near ties — iterations in which the runner-up term had a p-value within two-fold of the selected term — occurred in {int(ntt.sum())} of {len(nt)}
iterations of the full analysis; terms selected in a near tie were kept less often in the subsamples (median recurrence
{nt[ntt].recurrence_kept.median():.2f}, versus {nt[~ntt].recurrence_kept.median():.2f} for the other terms). Path dependence therefore mainly affects which
term represents a group of genes, not the significance or the gene composition of the network. Terms kept in fewer than half of the subsamples:</p>
<div class='tbl'><table><thead><tr><th>Collection</th><th>Term</th><th class='n'>q, full list</th><th class='n'>Kept in subsamples</th></tr></thead><tbody>{lowrows}</tbody></table></div>

<h2>C. Redundancy between collections</h2>
<p class="muted" style="font-size:13.5px">For each collection (rows), the number of its gene sets (10–600 genes after restriction to the background) that have a
Jaccard similarity of at least 0.9 with a gene set in another collection (columns).</p>
<div class='tbl'><table><thead><tr><th>Collection</th><th class='n'>Gene sets</th>{mhead}<th class='n'>Any other<br><span>J ≥ 0.9</span></th><th class='n'>Any other<br><span>identical</span></th></tr></thead><tbody>{mrows}</tbody></table></div>
<p class="notes">Deduplication: collections were ordered Hallmark, Reactome, KEGG Legacy, KEGG MEDICUS, BioCarta, PID, WikiPathways, GO BP, GO CC, GO MF, and a
gene set was removed if it had J ≥ 0.9 with a set in an earlier collection ({n_drop} sets removed: {E(drop_txt)}). Leave-one-out: each collection was
removed in turn. The HIV list was re-analyzed with 1,000 random lists per configuration (smallest attainable p = 0.001).</p>
<div class='tbl'><table><thead>{CHEAD}</thead><tbody>{crows}</tbody></table></div>

<h2>D. Annotation-matched random gene lists</h2>
<p class="muted" style="font-size:13.5px">Genes reported in screens and signatures tend to be well studied and therefore annotated to more gene sets than a random
background gene (mean number of the 10,439 analysed gene sets containing a gene: {drows[0]['ann_bg']:.1f} in the background). Random lists drawn uniformly
from the background could therefore underestimate the number of terms selected by chance. In this analysis each random list was matched to the input
gene by gene: for every input gene, a background gene was drawn from the same annotation bin (genes in no analysed gene set, plus nine quantile bins of
the annotation count of the annotated genes). The reduced-size lists for the q-value adjustment were random subsets of matched lists. Everything else
(10,000 lists, seed 2026, selection, q-values, term filter and network test) was as in the main analysis. Arrows: standard random lists → annotation-matched lists.</p>
<div class='tbl'><table><thead><tr><th>Gene list</th><th class='n'>Input annotation<br><span>mean gene sets per gene</span></th><th class='n'>Chance terms<br><span>per random list</span></th>
<th class='n'>Terms kept<br><span>q &lt; 0.10</span></th><th class='n'>Matched network<br><span>cluster genes / terms</span></th><th>Network p<br><span>cluster size</span></th>
<th>Network p<br><span>edges per gene</span></th></tr></thead><tbody>{DROWS}</tbody></table></div>
<p class="notes">The networks remain significant against annotation-matched lists, but random lists of well-annotated genes select more terms by chance,
so fewer terms pass the q-value filter. Annotation density therefore affects which individual terms are kept more than the significance of the network.
In the HIV network every community keeps at least one term:</p>
<div class='tbl'><table class='narrow'><thead><tr><th>Community</th><th class='n'>Terms kept<br><span>standard lists</span></th><th class='n'>Terms kept<br><span>matched lists</span></th></tr></thead><tbody>{CROWS}</tbody></table></div>
<p class="notes">The {len(dropped)} HIV terms not kept against matched lists had q-values between {dropped.q.min():.2f} and {dropped.q.max():.2f} in that analysis and
{dropped["q-value"].min():.3f}–{dropped["q-value"].max():.3f} in the main analysis, i.e. they were the weaker of the kept terms (every term with q &lt; {dropped["q-value"].min():.3f}
was kept in both analyses); they include all {n_unadj} terms whose q-value
without the signal adjustment is 0.10 or more (Supplementary Table 5). Individual terms near the threshold should therefore be interpreted with caution.</p>

<h2>Notes</h2><ul class="notes">
<li><b>Default analysis:</b> hypergeometric test (upper tail); p &lt; 0.01 for selection; at least 3 input genes
per term; gene sets of 10–600 genes after restriction to the background; at most 30 iterations per collection; the ten collections of Supplementary Table 2.</li>
<li><b>Random gene lists:</b> Part A, 2,000 lists of the input size and 200 per reduced size for the q-value adjustment (seed 2026), generated with the
setting of each row; Part B, 5,000 lists; Part C, 1,000 lists; Part D, 10,000 annotation-matched lists and 1,000 per reduced size. Empirical p = (r + 1)/(n + 1). The main analyses use 10,000 lists.</li>
<li><b>Terms kept:</b> terms with signal-adjusted empirical q &lt; 0.10, which form the network (an empirical estimate used to filter terms, not formal false discovery rate control). <b>Jaccard:</b> size of the intersection divided by the size of the union.</li>
</ul>
<footer>Code and data: iGEA repository, publication/csbj_2026/analyses — step10_extras/s10a_sensitivity.py (A), step06_stability/s06_stability.py (B),
step03_redundancy/s03_redundancy.py and s03b_sensitivity.py (C), step10_extras/s10c_annotation_matched.py (D). Builder: publication/csbj_2026/build/supplement/build_supp_table_08_robustness.py.</footer>
</main></body></html>"""
open(f"{OUT}/Supp Table 08 - Robustness analyses.html", "w").write(doc)
# ---------------------------------------------------------------- TXT
L = ["Supplementary Table 8. Robustness of iGEA results", "",
     "A. Analysis settings (one setting changed at a time; random lists generated with the same setting; 2,000 lists).", "",
     d.to_csv(sep="\t", index=False),
     "B. Input subsampling (100 random subsamples of 90% of the HIV list, 238 genes; 5,000 random lists).", ""]
L += [f"{a}\t{b}\t(full list: {c})" for a, b, c in B_ROWS]
L += ["", f"Full-list kept terms kept in >= 80% / >= 50% of subsamples: {int((k.recurrence_kept >= 0.8).sum())} / {int((k.recurrence_kept >= 0.5).sum())} of 56; "
      f"median recurrence {k.recurrence_kept.median():.2f}.",
      f"Near ties (runner-up p within two-fold): {int(ntt.sum())} of {len(nt)} iterations; median recurrence {nt[ntt].recurrence_kept.median():.2f} vs {nt[~ntt].recurrence_kept.median():.2f}.", "",
      "Terms kept in fewer than half of the subsamples:", low[["Library", "Term", "q_full", "recurrence_kept"]].to_csv(sep="\t", index=False),
      "C. Redundancy between collections (gene sets with Jaccard >= 0.9 to a set in another collection).", "", mat.to_csv(sep="\t", index=False),
      f"Deduplication removed {n_drop} sets ({drop_txt}). Re-analysis (1,000 random lists):", "",
      pd.concat([dedup.assign(Configuration="all 10 (deduplicated)"), loo])[["Configuration", "iGEA terms", "Kept (q < 0.10)", "Cluster genes", "Cluster terms",
       "Edges/gene", "p Cluster Size", "p Average Edges per Gene"]].to_csv(sep="\t", index=False),
      "D. Annotation-matched random gene lists (10,000 lists matched to the input gene by gene on annotation count; standard -> matched).", "",
      pd.DataFrame(drows).rename(columns={"ds": "Gene list", "sel": "Terms selected", "ann_in": "Input mean annotation", "ann_bg": "Background mean annotation",
        "ch_u": "Chance terms/list standard", "ch_m": "Chance terms/list matched", "kept_u": "Kept standard", "kept_m": "Kept matched",
        "p_cs_u": "p cluster size standard", "p_cs_m": "p cluster size matched", "p_epg_u": "p edges/gene standard", "p_epg_m": "p edges/gene matched",
        "genes_m": "Matched cluster genes", "terms_m": "Matched cluster terms"}).to_csv(sep="\t", index=False),
      "HIV terms kept per community (standard / matched):", cg.to_csv(sep="\t")]
open(f"{OUT}/Supp Table 08 - Robustness analyses.txt", "w").write("\n".join(L))
print("written")
