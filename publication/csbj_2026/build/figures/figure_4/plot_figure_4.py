"""Figure: iGEA networks of three independent HIV dependency-factor screens; shared themes in the same colour."""
import pathlib as _pl  # [publication] repository-relative paths
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys
sys.path.insert(0, str(_REPO / "code" / "ui"))
from network_plot import _layout, _hull
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
S9 = str(_PUB / "analyses" / "step09_end_to_end") + "/"
D = str(_PUB / "analyses" / "step11_three_screens") + "/"
THEME = {"Mediator": "#2CA02C", "Ubiquitin": "#8C564B", "RNA": "#9467BD", "Pore": "#FF7F0E", "Signalling": "#C51B7D",
         "own1": "#7F7F7F", "own2": "#A6A6A6", "own3": "#C7C7C7"}
THEME_NAME = {"Mediator": "Mediator / transcription coregulators", "Ubiquitin": "Ubiquitin–proteasome",
              "RNA": "RNA processing / RNA polymerase", "Pore": "Nuclear pore / export", "Signalling": "NF-κB / AKT signaling"}
PANELS = [
    (S9 + "HIV_dependency_factors.symbols_20261120_114347/", "A  Brass et al.", 265,
     {"C1": ("own1", "Golgi transport"), "C2": ("Signalling", "NF-κB/AKT signaling"), "C3": ("Mediator", "Mediator"),
      "C4": ("Pore", "nuclear pore"), "C5": ("RNA", "RNA exosome"), "C6": ("own2", "ER glycosylation, GPI"),
      "C7": ("Ubiquitin", "cullin-RING ligases")}),
    (D + "Konig_2008_HIV_dependency_factors.symbols_20261120_201945/", "B  König et al.", 292,
     {"C1": ("RNA", "pre-mRNA splicing"), "C2": ("Ubiquitin", "proteasome"), "C3": ("Mediator", "Mediator"),
      "C4": ("Pore", "nuclear pore"), "C5": ("RNA", "RNA polymerase II"),
      "C6": ("own1", "nucleases"), "C7": ("own2", "tyrosine phosphatases"), "C8": ("own3", "perinuclear ER")}),
    (D + "Zhou_2008_HIV_dependency_factors.symbols_20261120_202031/", "C  Zhou et al.", 311,
     {"C1": ("Signalling", "NF-κB/AKT signaling"), "C2": ("RNA", "RNA polymerase II"), "C3": ("Mediator", "Mediator"),
      "C4": ("own1", "mitochondrial complex I"), "C5": ("Ubiquitin", "K48 ubiquitination")}),
]
split = lambda s: [g.strip() for g in str(s).split(",") if g.strip()]


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


def draw(ax, d, title, n_in, names, seed=2026):
    x_fixed = None
    t = pd.read_csv(d + "combined_iterative_results.tsv", sep="\t"); t = t[t["In network"] == True]
    terms = [dict(library=r.Library, term=r.Term, genes=split(r.Genes)) for r in t.itertuples()]
    mt = pd.read_csv(d + "network_modules_terms.tsv", sep="\t")
    tc = {(r.Library, r.Term): r.Community for r in mt.itertuples()}
    mg = pd.read_csv(d + "network_modules_genes.tsv", sep="\t").set_index("Gene")
    gene_info = {g: r for g, r in zip(mg.index, mg.itertuples())}
    pos = _layout(terms, tc, gene_info, seed=seed)
    pos = pull_in(pos, terms)
    col = lambda c: THEME[names[c][0]]
    tcomm = {i: tc[(x["library"], x["term"])] for i, x in enumerate(terms)}
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    span = max(max(xs) - min(xs), max(ys) - min(ys))
    comms = sorted(set(tcomm.values()), key=lambda c: int(c[1:]))
    for c in comms:
        pts = [pos[("T", i)] for i, cc in tcomm.items() if cc == c] + \
              [pos[("G", g)] for g in mg.index if mg.loc[g, "Community"] == c and not mg.loc[g, "Bridge"] and ("G", g) in pos]
        poly = _hull(pts, 0.03 * span)
        ax.fill([p[0] for p in poly], [p[1] for p in poly], color=col(c), alpha=0.10, lw=0)
    edges = [(("T", i), ("G", g)) for i, x in enumerate(terms) for g in set(x["genes"])]
    for a, b in edges:
        ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#B0B0B0", lw=0.5, alpha=0.45, zorder=1)
    genes = sorted({g for x in terms for g in x["genes"]})
    ax.scatter([pos[("G", g)][0] for g in genes], [pos[("G", g)][1] for g in genes], s=3, zorder=2,
               c=[col(mg.loc[g, "Community"]) for g in genes],
               edgecolors=["black" if mg.loc[g, "Bridge"] else "none" for g in genes], linewidths=0.5)
    ax.scatter([pos[("T", i)][0] for i in tcomm], [pos[("T", i)][1] for i in tcomm], marker="s", zorder=3,
               s=[6 + 0.8 * len(x["genes"]) for x in terms], c=[col(tcomm[i]) for i in tcomm],
               edgecolors="white", linewidths=0.5, alpha=0.95)
    x0 = ax.get_position().x0; F = ax.figure
    F.text(x0, 0.955, f"({title[0].lower()})", fontsize=9, weight="bold", va="bottom", ha="left")
    F.text(x0 + 0.032, 0.955, title[3:], fontsize=8, va="bottom", ha="left")
    F.text(x0 + 0.032, 0.915, f"{n_in} input genes, {len(terms)} terms", fontsize=8, va="bottom", ha="left")
    ax.set_aspect("equal"); ax.set_anchor("N"); ax.axis("off")


plt.rcParams.update({"font.family": "Liberation Serif", "pdf.fonttype": 42, "ps.fonttype": 42})
fig = plt.figure(figsize=(7.48, 4.3))
for k, (d, title, n, names) in enumerate(PANELS):
    draw(fig.add_axes([0.01 + k * 0.333, 0.42, 0.32, 0.47]), d, title, n, names)
ax = fig.add_axes([0.01, 0.0, 0.98, 0.37])
# R2 (Dec 2026): FDR = Benjamini-Hochberg over the 131 cross-screen pairs, matched random communities (step11 s11h, K = 25);
# ranges over the best pair of each pair of screens (Supplementary Table 14). R1 showed uncorrected p (2,000 random groups).
rows = [("Mediator", "Mediator / transcription coregulators", ["C3", "C3", "C3"], "0.56–0.78", "< 0.001"),
        ("Ubiquitin", "Ubiquitin–proteasome", ["C7", "C2", "C5"], "0.08–0.16", "0.009–0.03"),
        ("RNA", "RNA processing / RNA polymerase", ["C5", "C1 C5", "C2"], "0.10–0.12", "0.05–0.13"),
        ("Pore", "Nuclear pore / export", ["C4", "C4", ""], "0.35", "0.05"),
        ("Signalling", "NF-κB / AKT signaling", ["C2", "", "C1"], "0.10", "0.90")
        ]
xs = [0.42, 0.53, 0.64]; jx, px = 0.79, 0.92
ax.set_xlim(0, 1); ax.set_ylim(len(rows) + 0.6, -1.0); ax.axis("off")
ax.text(0.0, -0.8, "(d)", fontsize=9, weight="bold", va="center"); ax.text(0.034, -0.8, "Themes shared between screens", fontsize=8, va="center")
for x, h in zip(xs + [jx, px], ["Brass", "König", "Zhou", "Jaccard", "FDR"]):
    ax.text(x, 0.0, h, ha="center", va="center", fontsize=7.5, weight="bold")
for i, (th, name, cs, j, p) in enumerate(rows, start=1):
    shared = not th.startswith("own")
    if i == 6:
        ax.plot([0, 1], [i - 0.55, i - 0.55], color="#999", lw=0.6)
        ax.text(0.79, i, "found in one screen only", ha="left", va="center", fontsize=7, color="#666", style="italic")

    if shared: ax.add_patch(plt.Rectangle((0, i - 0.45), 1, 0.9, color=THEME[th], alpha=0.10, lw=0))
    ax.text(0.01, i, name, fontsize=7.5, va="center", color="#111" if shared else "#444")
    for x, c in zip(xs, cs):
        parts = c.split()
        for k, lab in enumerate(parts):
            xx = x + (k - (len(parts) - 1) / 2) * 0.022
            ax.scatter([xx], [i], s=55, color=THEME[th], zorder=4, lw=0)
    ax.text(jx, i, j, ha="center", va="center", fontsize=7.5); ax.text(px, i, p, ha="center", va="center", fontsize=7.5)
fig.savefig("Figure 4 - Three HIV screens.pdf"); fig.savefig("Figure 4 - Three HIV screens.png", dpi=600)
