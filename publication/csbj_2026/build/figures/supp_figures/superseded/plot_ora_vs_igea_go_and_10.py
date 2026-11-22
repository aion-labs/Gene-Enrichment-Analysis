"""HIV networks, ORA vs iGEA, with 3 GO collections and with 10 collections; colours = themes (same across panels)."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import pandas as pd, networkx as nx, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B = str(_PUB) + "/analyses/step09_end_to_end/"
THEME = {"Golgi": "#1F77B4", "Signalling": "#D62728", "Mediator": "#2CA02C", "Pore": "#FF7F0E",
         "RNA": "#9467BD", "GPI": "#17BECF", "Cullin": "#8C564B"}
PANELS = [  # (folder, kind, title, community -> (theme, label))
    ("ora_go_fdr0.05/", "ora", "A. ORA, 3 GO collections",
     {"C1": ("RNA", "RNA processing, complex disassembly"), "C2": ("Pore", "nuclear pore, protein localization"),
      "C3": ("Golgi", "Golgi transport"), "C4": ("GPI", "glycolipid / GPI"), "C5": ("Mediator", "Mediator")}),
    ("hiv_go_only/HIV_dependency_factors.symbols_20261120_114431/", "igea", "B. iGEA, 3 GO collections",
     {"C1": ("RNA", "RNA processing, signalling"), "C2": ("Golgi", "Golgi transport"), "C3": ("Pore", "nuclear pore"),
      "C4": ("Mediator", "Mediator"), "C5": ("GPI", "ER glycosylation, GPI"), "C6": ("Cullin", "cullin-RING ligases")}),
    ("ora_all_fdr0.05/", "ora", "C. ORA, 10 collections",
     {"C1": ("Golgi", "Golgi, glycosylation, signalling, endocytosis"),
      "C2": ("Pore", "nuclear pore + viral, SUMO, mitosis, glucokinase terms")}),
    ("HIV_dependency_factors.symbols_20261120_114347/", "igea", "D. iGEA, 10 collections",
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


def packed_layout(G, seed):
    """Spring layout per connected component; components scaled by size and packed left to right, small ones stacked."""
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    pos, x0 = {}, 0.0
    big = comps[0]; H = G.subgraph(big)
    p = nx.spring_layout(H, seed=seed, k=0.9 / len(H) ** 0.5, iterations=300)
    pos.update({n: (v[0] * 1.0, v[1] * 1.0) for n, v in p.items()})
    x0 = 1.25; y0 = 0.9
    for c in comps[1:]:
        H = G.subgraph(c); r = 0.12 + 0.5 * (len(c) / len(big)) ** 0.5
        p = nx.spring_layout(H, seed=seed, k=0.9 / len(H) ** 0.5, iterations=200) if len(H) > 2 else {n: (0, 0) for n in H}
        pos.update({n: (x0 + v[0] * r, y0 - r + v[1] * r) for n, v in p.items()})
        y0 -= 2.2 * r + 0.1
    return pos


def draw(ax, folder, kind, title, names, seed=7):
    terms, tc, mg = load(folder, kind)
    G = nx.Graph()
    for (lib, t), genes in terms.items():
        G.add_node(("T", t), kind="term", comm=tc[(lib, t)], size=len(genes))
        for g in genes:
            G.add_node(("G", g), kind="gene", comm=mg.loc[g, "Community"], bridge=bool(mg.loc[g, "Bridge"]))
            G.add_edge(("T", t), ("G", g))
    col = lambda c: THEME[names[c][0]]
    pos = packed_layout(G, seed)
    nx.draw_networkx_edges(G, pos, ax=ax, width=0.25, alpha=0.25, edge_color="#777")
    genes = [n for n, d in G.nodes(data=True) if d["kind"] == "gene"]
    trms = [n for n, d in G.nodes(data=True) if d["kind"] == "term"]
    nx.draw_networkx_nodes(G, pos, nodelist=genes, ax=ax, node_size=12, node_color=[col(G.nodes[n]["comm"]) for n in genes],
                           edgecolors=["black" if G.nodes[n]["bridge"] else "none" for n in genes], linewidths=0.6)
    nx.draw_networkx_nodes(G, pos, nodelist=trms, ax=ax, node_shape="s", node_size=[18 + 3 * G.nodes[n]["size"] for n in trms],
                           node_color=[col(G.nodes[n]["comm"]) for n in trms], edgecolors="white", linewidths=0.5, alpha=0.9)
    deg = mg["Degree"]
    comms = sorted({d["comm"] for _, d in G.nodes(data=True)}, key=lambda c: int(c[1:]))
    ax.set_title(f"{title}\n{len(trms)} terms · {len(genes)} genes · {G.number_of_edges()} edges · "
                 f"{len(comms)} communities · max gene degree {deg.max()}", fontsize=10)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=col(c), label=f"{c} {names[c][1]}") for c in comms],
              loc="upper left", bbox_to_anchor=(0, -0.0), fontsize=7.5, frameon=False, ncol=2, handletextpad=0.2, columnspacing=1)
    ax.axis("off")


fig, axs = plt.subplots(2, 2, figsize=(15, 15.5))
for ax, (folder, kind, title, names) in zip(axs.flat, PANELS):
    draw(ax, folder, kind, title, names)
fig.text(0.5, 0.008, "HIV dependency factors (265 genes). ORA: terms with BH FDR < 0.05 (within each collection), each linked to all "
         "its input genes. iGEA: terms with empirical q < 0.10, each gene assigned to at most one term per collection.\n"
         "Squares = terms (size = number of genes), circles = genes, black outline = bridge gene (terms in more than one "
         "community). Communities: Leiden on the term network (same procedure for all panels); colours = shared themes.",
         ha="center", fontsize=8.5, color="#333")
plt.tight_layout(rect=(0, 0.035, 1, 1), h_pad=4)
plt.savefig("out/ORA_vs_iGEA_GO_and_10_collections.png", dpi=200)
plt.savefig("out/ORA_vs_iGEA_GO_and_10_collections.pdf")
