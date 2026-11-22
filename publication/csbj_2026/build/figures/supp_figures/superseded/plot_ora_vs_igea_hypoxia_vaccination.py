"""Supplementary figure: hypoxia and vaccination networks with 10 collections, ORA vs iGEA (app layout)."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys
sys.path.insert(0, str(_REPO) + "/code/ui")
from network_plot import _layout, _hull
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
S = str(_PUB) + "/build/figures/supp_figure_6_data/"
PAL = ["#1F77B4", "#FF7F0E", "#2CA02C", "#9467BD", "#17BECF", "#8C564B", "#E377C2", "#7F7F7F"]
split = lambda s: [g.strip() for g in str(s).split(",") if g.strip()]
def nice(t):
    t = t.replace("_", " ")
    for a, b in [("GOBP ", "GO BP: "), ("GOCC ", "GO CC: "), ("GOMF ", "GO MF: "), ("REACTOME ", "Reactome: "), ("HALLMARK ", "Hallmark: ")]:
        if t.startswith(a): t = b + t[len(a):]
    t = t.replace("GOBP: ", "GO BP: ").replace("GOCC: ", "GO CC: ").replace("GOMF: ", "GO MF: ").replace("REACTOME: ", "Reactome: ").replace("HALLMARK: ", "Hallmark: ")
    head, _, rest = t.partition(": ")
    return head + ": " + rest.lower() if rest else t.lower()
def draw(ax, d, kind, title):
    tt = pd.read_csv(S + d + ("/ora_terms.tsv" if kind == "ora" else "/igea_terms.tsv"), sep="\t")
    terms = [dict(library=r.Collection, term=r.Term, genes=split(r.Genes)) for r in tt.itertuples()]
    mt = pd.read_csv(S + d + "/network_modules_terms.tsv", sep="\t"); mg = pd.read_csv(S + d + "/network_modules_genes.tsv", sep="\t").set_index("Gene")
    mc = pd.read_csv(S + d + "/network_modules_communities.tsv", sep="\t")
    tc = {(r.Library, r.Term): r.Community for r in mt.itertuples()}
    gene_info = {g: r for g, r in zip(mg.index, mg.itertuples())}
    pos = _layout(terms, tc, gene_info, seed=2026)
    comms = list(mc.Community); col = {c: PAL[i % len(PAL)] for i, c in enumerate(comms)}
    tcomm = {i: tc[(t["library"], t["term"])] for i, t in enumerate(terms)}
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]; span = max(max(xs) - min(xs), max(ys) - min(ys))
    for c in comms:
        pts = [pos[("T", i)] for i, cc in tcomm.items() if cc == c] + [pos[("G", g)] for g in mg.index if mg.loc[g, "Community"] == c and not mg.loc[g, "Bridge"] and ("G", g) in pos]
        poly = _hull(pts, 0.03 * span); ax.fill([p[0] for p in poly], [p[1] for p in poly], color=col[c], alpha=0.08, lw=0)
    edges = [(("T", i), ("G", g)) for i, t in enumerate(terms) for g in set(t["genes"])]
    for a, b in edges:
        ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#777", lw=0.2, alpha=0.25, zorder=1)
    genes = sorted({g for t in terms for g in t["genes"]})
    ax.scatter([pos[("G", g)][0] for g in genes], [pos[("G", g)][1] for g in genes], s=10, zorder=2, c=[col[mg.loc[g, "Community"]] for g in genes],
               edgecolors=["black" if mg.loc[g, "Bridge"] else "none" for g in genes], linewidths=0.5)
    ax.scatter([pos[("T", i)][0] for i in tcomm], [pos[("T", i)][1] for i in tcomm], marker="s", zorder=3, s=[16 + 2 * len(t["genes"]) for t in terms],
               c=[col[tcomm[i]] for i in tcomm], edgecolors="white", linewidths=0.4, alpha=0.95)
    deg = mg["Degree"]; s = dict(l.rstrip("\n").split("\t") for l in open(S + d + "/network_modules_summary.tsv"))
    ax.set_title(f"{title}\n{len(terms)} terms · {len(genes)} genes · {len(edges)} edges · {len(comms)} communities "
                 f"(modularity {float(s['modularity']):.2f}) · max gene degree {deg.max()}", fontsize=9.5)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=col[r.Community], label=f"{r.Community} ({r.Terms} terms) {nice(r._5)[:48]}") for r in mc.itertuples()],
              loc="upper left", bbox_to_anchor=(0, 0), fontsize=7, frameon=False, ncol=1, handletextpad=0.2)
    ax.set_aspect("equal"); ax.axis("off")
fig, axs = plt.subplots(2, 2, figsize=(15, 16))
for ax, (d, k, t) in zip(axs.flat, [("Hypoxia_ORA_all", "ora", "A. Hypoxia (75 genes), ORA, 10 collections"),
                                     ("Hypoxia_iGEA_all", "igea", "B. Hypoxia (75 genes), iGEA, 10 collections"),
                                     ("Vaccination_ORA_all", "ora", "C. Vaccination (426 genes), ORA, 10 collections"),
                                     ("Vaccination_iGEA_all", "igea", "D. Vaccination (426 genes), iGEA, 10 collections")]):
    draw(ax, d, k, t)
fig.text(0.5, 0.006, "ORA: terms with BH FDR < 0.05 (within each collection; gene sets of 10–600 genes, as for iGEA), each linked to all its input genes. "
         "iGEA: terms with empirical q < 0.10, each gene assigned to at most one term per collection.\nSquares = terms (size = number of genes), circles = genes, "
         "black outline = bridge gene. Communities: Leiden on the term network (same procedure for all panels), shaded; legend = community and its lead term. Layout as in the iGEA app.",
         ha="center", fontsize=8, color="#333")
plt.tight_layout(rect=(0, 0.03, 1, 1), h_pad=7)
plt.savefig("out/ORA_vs_iGEA_hypoxia_vaccination_10_collections.png", dpi=200)
plt.savefig("out/ORA_vs_iGEA_hypoxia_vaccination_10_collections.pdf")
