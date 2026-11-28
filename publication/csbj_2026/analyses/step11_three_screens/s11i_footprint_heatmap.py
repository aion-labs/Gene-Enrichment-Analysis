"""Footprint Jaccard between all communities of the three HIV screens (iGEA networks), clustered.
Footprint = union of the full gene sets (background-restricted, 10-600 genes) of a community's terms.
Rows/columns: 'screen:community · lead term'; label colors = screen. Average-linkage clustering on 1 - Jaccard.
Output: s11i_footprint_jaccard_heatmap.pdf/.png, s11i_footprint_jaccard_matrix.tsv"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns, matplotlib.lines
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
H = os.path.dirname(os.path.abspath(__file__))
ns = {"__file__": os.path.join(H, "s11b_compare.py"), "__name__": "x"}
exec(open(os.path.join(H, "s11b_compare.py")).read().split("rows = []")[0], ns)
members, data = ns["members"], ns["data"]
COL = {"Brass": "#1f77b4", "König": "#ff7f0e", "Zhou": "#2ca02c"}
def lead(s):
    s = s.split(": ", 1)[-1].lower()
    for a, b in (("rna", "RNA"), ("dna", "DNA"), ("gtpase", "GTPase"), ("k48", "K48"), ("mrna", "mRNA"), ("pre mrna", "pre-mRNA")):
        s = s.replace(a, b)
    return s if len(s) <= 38 else s[:36] + "…"
# short community names as in Supplementary Table 14 and Figure 4
SHORT = {("Brass", "C1"): "Golgi transport", ("Brass", "C2"): "EGFR/AKT/NF-κB signaling", ("Brass", "C3"): "Mediator",
         ("Brass", "C4"): "nuclear pore, export", ("Brass", "C5"): "RNA exosome, RNA processing", ("Brass", "C6"): "ER glycosylation, GPI anchor",
         ("Brass", "C7"): "cullin-RING ligases", ("König", "C1"): "mRNA splicing", ("König", "C2"): "proteasome",
         ("König", "C3"): "Mediator, transcription coregulators", ("König", "C4"): "nuclear pore, import", ("König", "C5"): "RNA polymerase II, DNA repair",
         ("König", "C6"): "nucleases", ("König", "C7"): "tyrosine phosphatases", ("König", "C8"): "perinuclear ER",
         ("Zhou", "C1"): "NF-κB/AKT signaling, autophagy", ("Zhou", "C2"): "RNA polymerase II, entry receptors", ("Zhou", "C3"): "Mediator",
         ("Zhou", "C4"): "mitochondrial complex I", ("Zhou", "C5"): "K48 ubiquitination"}
names, fps = [], []
for s in data:
    for c in sorted(data[s]["comm"], key=lambda x: int(x[1:])):
        names.append((s, c)); fps.append(set().union(*[members[t] for t in data[s]["comm"][c]]))
n = len(names)
J = np.array([[len(a & b) / len(a | b) for b in fps] for a in fps])
labels = [f"{s}:{c} · {SHORT[(s, c)]}" for s, c in names]
pd.DataFrame(J, index=labels, columns=labels).round(3).to_csv(os.path.join(H, "s11i_footprint_jaccard_matrix.tsv"), sep="\t")
Z = linkage(squareform(1 - J, checks=False), method="average")
show = pd.DataFrame(J, index=labels, columns=labels)
# themes: cut the dendrogram at 1 - Jaccard = 0.92 (average Jaccard >= 0.08 between merged groups); clusters with
# communities from at least two screens are themes, the rest are screen-specific
from scipy.cluster.hierarchy import fcluster
CUT = 0.92
cl = fcluster(Z, t=CUT, criterion="distance")
THEME_NAMES = {frozenset({"Brass:C3", "König:C3", "Zhou:C3"}): "Mediator / transcription coregulators",
               frozenset({"Brass:C7", "König:C2", "Zhou:C5"}): "ubiquitin–proteasome",
               frozenset({"Brass:C1", "Brass:C4", "König:C4"}): "nuclear pore and transport",
               frozenset({"Brass:C2", "Zhou:C1"}): "NF-κB / AKT signaling",
               frozenset({"Brass:C5", "König:C1", "Zhou:C2"}): "RNA processing and RNA polymerase"}
THEME_COL = ["#d62728", "#9467bd", "#8c564b", "#e377c2", "#bcbd22"]   # no blue, orange or green (screen colors)
themes = {}
for k in sorted(set(cl)):
    idx = [i for i in range(n) if cl[i] == k]
    if len({names[i][0] for i in idx}) >= 2:
        key = frozenset(f"{names[i][0]}:{names[i][1]}" for i in idx)
        themes[k] = (THEME_NAMES.get(key, " + ".join(sorted(key))), idx)
order_themes = sorted(themes, key=lambda k: -max(J[a, b] for a in themes[k][1] for b in themes[k][1] if a != b))
tcol = {k: THEME_COL[j % len(THEME_COL)] for j, k in enumerate(order_themes)}
theme_of = {i: k for k, (_, idx) in themes.items() for i in idx}
side = pd.DataFrame({"screen": [COL[s] for s, _ in names],
                     "theme": [tcol[theme_of[i]] if i in theme_of else "#ffffff" for i in range(n)]}, index=labels)
mask = np.eye(n, dtype=bool)
# adjusted significance of each pair (matched null, K = 25, BH over the 131 iGEA pairs; s11h): * FDR < 0.05, † FDR < 0.10
mt = pd.read_csv(os.path.join(H, "s11h_igea_pairs_matched.tsv"), sep="\t")
fdr = {}
for r in mt.itertuples():
    fdr[(r.screen_1, r.community_1, r.screen_2, r.community_2)] = fdr[(r.screen_2, r.community_2, r.screen_1, r.community_1)] = r.matched_BH_FDR
def mark(i, j):
    q = fdr.get((*names[i], *names[j]))
    return "" if q is None else "*" if q < 0.05 else "†" if q < 0.10 else ""
ann = np.array([[f"{J[i, j]:.2f}{mark(i, j)}" for j in range(n)] for i in range(n)])
g = sns.clustermap(show, row_linkage=Z, col_linkage=Z, cmap="YlOrRd", vmin=0, vmax=0.4, mask=mask,
                   row_colors=side, col_colors=side, figsize=(16, 15), dendrogram_ratio=(0.10, 0.10),
                   colors_ratio=0.025, cbar_pos=(0.02, 0.86, 0.02, 0.10), cbar_kws={"label": "footprint Jaccard (≥ 0.4 saturated)"},
                   linewidths=0.5, linecolor="white", annot=ann, fmt="", annot_kws={"size": 6.5})
order = g.dendrogram_row.reordered_ind
pos = {o: p for p, o in enumerate(order)}
for t in g.ax_heatmap.get_xticklabels() + g.ax_heatmap.get_yticklabels():
    s_ = t.get_text().split(":")[0]; t.set_color(COL[s_]); t.set_fontweight("bold"); t.set_fontsize(8.5)
for txt in g.ax_heatmap.texts:
    try: v = float(txt.get_text().rstrip("*†"))
    except ValueError: continue
    txt.set_color("white" if v > 0.3 else ("#222222" if v >= 0.05 else "#c8c8c8"))
# outline each theme block on the heatmap and number it
from matplotlib.patches import Rectangle, Patch
for j, k in enumerate(order_themes, start=1):
    ps = sorted(pos[i] for i in themes[k][1]); a, b = ps[0], ps[-1] + 1
    assert ps == list(range(a, b))                      # dendrogram cut -> contiguous block
    g.ax_heatmap.add_patch(Rectangle((a, a), b - a, b - a, fill=False, ec=tcol[k], lw=3, clip_on=False))
    g.ax_heatmap.text(b + 0.1, a + 0.05, str(j), color=tcol[k], fontsize=13, fontweight="bold", va="top", ha="left")
g.ax_heatmap.set_xlabel(""); g.ax_heatmap.set_ylabel("")
for ax_ in (g.ax_row_colors, g.ax_col_colors):
    ax_.tick_params(labelsize=8)
leg1 = [Patch(color=v, label=k) for k, v in COL.items()]
l1 = g.ax_col_dendrogram.legend(handles=leg1, loc="lower left", ncol=1, frameon=False, fontsize=9.5,
                                bbox_to_anchor=(0.0, 1.0), title="Screens (labels, first bar)", title_fontsize=9.5)
g.ax_col_dendrogram.add_artist(l1)
g.ax_col_dendrogram.text(0.32, 1.06, "Themes (outlined, numbered blocks; second bar; named on the right):\n"
                         "clusters containing communities from at least two screens", transform=g.ax_col_dendrogram.transAxes,
                         fontsize=9.5, va="bottom")
g.fig.suptitle("Footprint Jaccard between the communities of the three HIV screens (iGEA, ten collections)\n"
               f"footprint = union of the full gene sets of a community's terms; average-linkage clustering on 1 − Jaccard, cut at {CUT} "
               "(average Jaccard ≥ 0.08); themes = clusters with communities from ≥ 2 screens; diagonal omitted\n"
               "* FDR < 0.05, † FDR < 0.10 vs random communities of matched size, Benjamini–Hochberg over the 131 cross-screen pairs",
               x=0.55, y=1.12, fontsize=11, fontweight="bold")
# theme brackets and names to the right of the row labels, spanning exactly the rows of each theme
import matplotlib.transforms as mtrans
g.fig.canvas.draw()
r = g.fig.canvas.get_renderer()
x_right = max(t.get_window_extent(r).x1 for t in g.ax_heatmap.get_yticklabels())
x_fig = g.fig.transFigure.inverted().transform((x_right, 0))[0] + 0.008
tr = mtrans.blended_transform_factory(g.fig.transFigure, g.ax_heatmap.transData)
for j, k in enumerate(order_themes, start=1):
    ps = sorted(pos[i] for i in themes[k][1]); a, b = ps[0] + 0.12, ps[-1] + 0.88
    c = tcol[k]
    g.fig.lines.append(matplotlib.lines.Line2D([x_fig, x_fig], [a, b], transform=tr, color=c, lw=4, solid_capstyle="butt"))
    for y in (a, b):
        g.fig.lines.append(matplotlib.lines.Line2D([x_fig - 0.006, x_fig], [y, y], transform=tr, color=c, lw=2))
    g.fig.text(x_fig + 0.008, (a + b) / 2, f"{j}  {themes[k][0]}", transform=tr, color=c, fontsize=11,
               fontweight="bold", va="center", ha="left")
out = os.path.join(H, "s11i_footprint_jaccard_heatmap")
g.savefig(out + ".pdf", bbox_inches="tight"); g.savefig(out + ".png", dpi=200, bbox_inches="tight")
pd.DataFrame([dict(theme=j, name=themes[k][0], communities="; ".join(labels[i] for i in themes[k][1]))
              for j, k in enumerate(order_themes, start=1)]).to_csv(os.path.join(H, "s11i_themes.tsv"), sep="\t", index=False)
print("saved", out, n, "communities;", len(themes), "themes")
