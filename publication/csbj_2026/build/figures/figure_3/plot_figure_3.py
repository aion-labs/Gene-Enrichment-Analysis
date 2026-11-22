"""HIV networks, ORA vs iGEA, with 3 GO collections and with 10 collections; colours = themes (same across panels)."""
import pathlib as _pl  # [publication] repository-relative paths
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys
sys.path.insert(0, str(_REPO / "code" / "ui"))
from network_plot import _layout, _hull   # the app's layout (community-aware) and community shading
import pandas as pd, networkx as nx, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B = str(_PUB / "analyses" / "step09_end_to_end") + "/"
THEME = {"Golgi": "#1F77B4", "Signalling": "#C51B7D", "Mediator": "#2CA02C", "Pore": "#FF7F0E",
         "RNA": "#9467BD", "GPI": "#17BECF", "Cullin": "#8C564B"}
PANELS = [  # (folder, kind, title, community -> (theme, label))
    ("ora_go_fdr0.05/", "ora", "A  ORA, 3 GO collections",
     {"C1": ("RNA", "RNA processing, complex disassembly"), "C2": ("Pore", "nuclear pore, protein localization"),
      "C3": ("Golgi", "Golgi transport"), "C4": ("GPI", "glycolipid / GPI"), "C5": ("Mediator", "Mediator")}),
    ("hiv_go_only/HIV_dependency_factors.symbols_20261120_114431/", "igea", "B  iGEA, 3 GO collections",
     {"C1": ("RNA", "RNA processing, signalling"), "C2": ("Golgi", "Golgi transport"), "C3": ("Pore", "nuclear pore"),
      "C4": ("Mediator", "Mediator"), "C5": ("GPI", "ER glycosylation, GPI"), "C6": ("Cullin", "cullin-RING ligases")}),
    ("ora_all_fdr0.05/", "ora", "C  ORA, 10 collections",
     {"C1": ("Golgi", "Golgi, glycosylation, signalling, endocytosis"),
      "C2": ("Pore", "nuclear pore + viral, SUMO, mitosis, glucokinase terms")}),
    ("HIV_dependency_factors.symbols_20261120_114347/", "igea", "D  iGEA, 10 collections",
     {"C1": ("Golgi", "Golgi transport"), "C2": ("Signalling", "EGFR/AKT/NF-κB, endocytosis"), "C3": ("Mediator", "Mediator"),
      "C4": ("Pore", "nuclear pore"), "C5": ("RNA", "RNA exosome"), "C6": ("GPI", "glycosylation, GPI"),
      "C7": ("Cullin", "cullin-RING ligases")}),
]
split = lambda s: [g.strip() for g in str(s).split(",") if g.strip()]


def load(folder, kind):
    d = B + folder
    if kind == "ora":
        t = pd.read_csv(d + "ora_terms.tsv", sep="\t")
        terms = {(r.Collection, r.Term): split(r.Genes) for r in t.itertuples()}
    else:
        t = pd.read_csv(d + "combined_iterative_results.tsv", sep="\t"); t = t[t["In network"] == True]
        terms = {(r.Library, r.Term): split(r.Genes) for r in t.itertuples()}
    mt = pd.read_csv(d + "network_modules_terms.tsv", sep="\t")
    mg = pd.read_csv(d + "network_modules_genes.tsv", sep="\t").set_index("Gene")
    return terms, {(r.Library, r.Term): r.Community for r in mt.itertuples()}, mg


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


TITLE_X = {"A": 0.01, "B": 0.51, "C": 0.01, "D": 0.51}


def draw(ax, folder, kind, title, names, seed=2026):
    terms_d, tc, mg = load(folder, kind)
    terms = [dict(library=l, term=t, genes=g) for (l, t), g in terms_d.items()]
    gene_info = {g: r for g, r in zip(mg.index, mg.itertuples())}
    pos = pull_in(_layout(terms, tc, gene_info, seed=seed), terms)
    col = lambda c: THEME[names[c][0]]
    tcomm = {i: tc[(t["library"], t["term"])] for i, t in enumerate(terms)}
    # community shading (as in the app): hull around each community's terms and non-bridge genes
    span = max(max(p[0] for p in pos.values()) - min(p[0] for p in pos.values()),
               max(p[1] for p in pos.values()) - min(p[1] for p in pos.values()))
    comms = sorted(set(tcomm.values()), key=lambda c: int(c[1:]))
    for c in comms:
        pts = [pos[("T", i)] for i, cc in tcomm.items() if cc == c] + \
              [pos[("G", g)] for g in mg.index if mg.loc[g, "Community"] == c and not mg.loc[g, "Bridge"] and ("G", g) in pos]
        poly = _hull(pts, 0.03 * span)
        ax.fill([p[0] for p in poly], [p[1] for p in poly], color=col(c), alpha=0.10, lw=0)
    edges = [(("T", i), ("G", g)) for i, t in enumerate(terms) for g in set(t["genes"])]
    for a, b in edges:
        ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#B0B0B0", lw=0.5, alpha=0.45, zorder=1)
    genes = sorted({g for t in terms for g in t["genes"]})
    ax.scatter([pos[("G", g)][0] for g in genes], [pos[("G", g)][1] for g in genes], s=3, zorder=2,
               c=[col(mg.loc[g, "Community"]) for g in genes],
               edgecolors=["black" if mg.loc[g, "Bridge"] else "none" for g in genes], linewidths=0.5)
    ax.scatter([pos[("T", i)][0] for i in tcomm], [pos[("T", i)][1] for i in tcomm], marker="s", zorder=3,
               s=[5 + 0.7 * len(t["genes"]) for t in terms], c=[col(tcomm[i]) for i in tcomm],
               edgecolors="white", linewidths=0.5, alpha=0.95)
    bb = ax.get_position()
    x0 = TITLE_X[title[0]]
    y2 = bb.y1 + 0.008; y1 = y2 + 0.027
    ax.figure.text(x0, y1, f"({title[0].lower()})", fontsize=9, fontweight="bold", va="bottom", ha="left")
    ax.figure.text(x0 + 0.03, y1, title[3:], fontsize=8, va="bottom", ha="left")
    ax.figure.text(x0 + 0.03, y2, f"{len(terms)} terms, {len(genes)} genes, {len(comms)} communities, "
                   f"{len(edges) / len(genes):.1f} edges per gene", fontsize=8, va="bottom", ha="left")
    ax.set_aspect("equal"); ax.set_anchor("C"); ax.axis("off")


plt.rcParams.update({"font.family": "Liberation Serif", "pdf.fonttype": 42, "ps.fonttype": 42})
fig = plt.figure(figsize=(7.48, 6.85))
boxes = [[0.01, 0.545, 0.48, 0.385], [0.51, 0.545, 0.48, 0.385], [0.01, 0.075, 0.48, 0.385], [0.51, 0.075, 0.48, 0.385]]
for bx, (folder, kind, title, names) in zip(boxes, PANELS):
    draw(fig.add_axes(bx), folder, kind, title, names)
from matplotlib.lines import Line2D as _L
for xs, ys in [((0.5, 0.5), (0.07, 0.995)), ((0.005, 0.995), (0.532, 0.532))]:
    fig.add_artist(_L(xs, ys, transform=fig.transFigure, color="#C8C8C8", lw=0.5))
KEY = [("Golgi", "Golgi transport"), ("Signalling", "Signaling (EGFR/AKT/NF-κB)"), ("Mediator", "Mediator"),
       ("Pore", "Nuclear pore"), ("RNA", "RNA processing"), ("GPI", "ER glycosylation, GPI"), ("Cullin", "Cullin-RING ligases")]
fig.legend(handles=[Line2D([], [], marker="s", ls="", color=THEME[k], markersize=6, label=l) for k, l in KEY],
           loc="lower center", ncol=4, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=1.5)
fig.savefig("Figure 3 - ORA vs iGEA networks.pdf"); fig.savefig("Figure 3 - ORA vs iGEA networks.png", dpi=600)
