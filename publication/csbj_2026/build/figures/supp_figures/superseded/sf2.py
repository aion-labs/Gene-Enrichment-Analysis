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
def pull_in(pos, terms):
    """Move small disconnected components next to the main one (the app spreads them far apart)."""
    import networkx as nx, numpy as np
    G = nx.Graph(); G.add_nodes_from(pos)
    G.add_edges_from((("T", i), ("G", g)) for i, x in enumerate(terms) for g in x["genes"])
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    P = {n: np.array(v, float) for n, v in pos.items()}
    main = np.array([P[n] for n in comps[0]]); c0 = main.mean(0)
    R = np.max(np.linalg.norm(main - c0, axis=1))
    for k, comp in enumerate(comps[1:]):
        pts = np.array([P[n] for n in comp]); cc = pts.mean(0); r = np.max(np.linalg.norm(pts - cc, axis=1)) if len(pts) > 1 else 0
        d = cc - c0; d = d / (np.linalg.norm(d) or 1)
        new = c0 + d * (R + r + 0.12 * R)
        for n in comp: P[n] = P[n] - cc + new
    return {n: tuple(v) for n, v in P.items()}


def draw(ax, d, kind, title):
    tt = pd.read_csv(S + d + ("/ora_terms.tsv" if kind == "ora" else "/igea_terms.tsv"), sep="\t")
    terms = [dict(library=r.Collection, term=r.Term, genes=split(r.Genes)) for r in tt.itertuples()]
    mt = pd.read_csv(S + d + "/network_modules_terms.tsv", sep="\t"); mg = pd.read_csv(S + d + "/network_modules_genes.tsv", sep="\t").set_index("Gene")
    mc = pd.read_csv(S + d + "/network_modules_communities.tsv", sep="\t")
    tc = {(r.Library, r.Term): r.Community for r in mt.itertuples()}
    gene_info = {g: r for g, r in zip(mg.index, mg.itertuples())}
    pos = pull_in(_layout(terms, tc, gene_info, seed=2026), terms)
    comms = list(mc.Community); col = {c: PAL[i % len(PAL)] for i, c in enumerate(comms)}
    tcomm = {i: tc[(t["library"], t["term"])] for i, t in enumerate(terms)}
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]; span = max(max(xs) - min(xs), max(ys) - min(ys))
    for c in comms:
        pts = [pos[("T", i)] for i, cc in tcomm.items() if cc == c] + [pos[("G", g)] for g in mg.index if mg.loc[g, "Community"] == c and not mg.loc[g, "Bridge"] and ("G", g) in pos]
        poly = _hull(pts, 0.03 * span); ax.fill([p[0] for p in poly], [p[1] for p in poly], color=col[c], alpha=0.10, lw=0)
    edges = [(("T", i), ("G", g)) for i, t in enumerate(terms) for g in set(t["genes"])]
    for a, b in edges:
        ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#777", lw=0.15, alpha=0.25, zorder=1)
    genes = sorted({g for t in terms for g in t["genes"]})
    ax.scatter([pos[("G", g)][0] for g in genes], [pos[("G", g)][1] for g in genes], s=3, zorder=2, c=[col[mg.loc[g, "Community"]] for g in genes],
               edgecolors=["black" if mg.loc[g, "Bridge"] else "none" for g in genes], linewidths=0.35)
    ax.scatter([pos[("T", i)][0] for i in tcomm], [pos[("T", i)][1] for i in tcomm], marker="s", zorder=3, s=[5 + 0.6 * len(t["genes"]) for t in terms],
               c=[col[tcomm[i]] for i in tcomm], edgecolors="white", linewidths=0.4, alpha=0.95)
    deg = mg["Degree"]; s = dict(l.rstrip("\n").split("\t") for l in open(S + d + "/network_modules_summary.tsv"))
    bb = ax.get_position()
    ax.figure.text(bb.x0, bb.y1 + 0.008, f"{title}\n{len(terms)} terms, {len(genes)} genes, {len(comms)} communities, "
                   f"{len(edges) / len(genes):.1f} edges per gene", fontsize=8, va="bottom", ha="left")
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=col[r.Community], label=f"{nice(r._5)[:52]} ({r.Terms})") for r in mc.itertuples()],
              loc="upper left", bbox_to_anchor=(0, -0.01), fontsize=7, frameon=False, ncol=1, handletextpad=0.2, labelspacing=0.2, borderaxespad=0)
    ax.set_aspect("equal"); ax.set_anchor("S"); ax.axis("off")
plt.rcParams.update({"font.family": "DejaVu Sans"})
fig = plt.figure(figsize=(7.48, 8.2))
pos_ = [[0.01, 0.635, 0.47, 0.30], [0.52, 0.635, 0.47, 0.30], [0.01, 0.15, 0.47, 0.30], [0.52, 0.15, 0.47, 0.30]]
for p_, (d, k, t) in zip(pos_, [("Hypoxia_ORA_all", "ora", "A  Hypoxia (75 genes), ORA"), ("Hypoxia_iGEA_all", "igea", "B  Hypoxia, iGEA"),
                                ("Vaccination_ORA_all", "ora", "C  Vaccination (426 genes), ORA"), ("Vaccination_iGEA_all", "igea", "D  Vaccination, iGEA")]):
    draw(fig.add_axes(p_), d, k, t)
fig.savefig("out/SF2.pdf"); fig.savefig("out/SF2.png", dpi=300)
