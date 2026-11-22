import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import pandas as pd, networkx as nx, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B = str(_PUB) + "/analyses/step09_end_to_end/"
PAL = ["#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD", "#17BECF", "#8C564B"]
split = lambda s: [g.strip() for g in str(s).split(",") if g.strip()]

def load(kind):
    if kind == "ora":
        d = B + "ora_all_fdr0.05/"; t = pd.read_csv(d + "ora_terms.tsv", sep="\t")
        terms = {(r.Collection, r.Term): split(r.Genes) for r in t.itertuples()}
    else:
        d = B + "HIV_dependency_factors.symbols_20261120_114347/"
        t = pd.read_csv(d + "combined_iterative_results.tsv", sep="\t"); t = t[t["In network"] == True]
        terms = {(r.Library, r.Term): split(r.Genes) for r in t.itertuples()}
    mt = pd.read_csv(d + "network_modules_terms.tsv", sep="\t")
    mg = pd.read_csv(d + "network_modules_genes.tsv", sep="\t")
    tc = {(r.Library, r.Term): r.Community for r in mt.itertuples()}
    return terms, tc, mg.set_index("Gene")

NAMES = {"ora": {"C1": "Golgi, glycosylation, signalling, endocytosis (33 terms)",
                 "C2": "nuclear pore, RNA export, viral infection, SUMOylation (54 terms)"},
         "igea": {"C1": "Golgi transport", "C2": "EGFR/AKT/NF-κB, endocytosis", "C3": "Mediator",
                  "C4": "nuclear pore", "C5": "RNA exosome", "C6": "glycosylation, GPI", "C7": "cullin-RING ligases"}}


def draw(ax, kind, title, seed=7):
    terms, tc, mg = load(kind)
    G = nx.Graph()
    for (lib, t), genes in terms.items():
        G.add_node(("T", t), kind="term", comm=tc[(lib, t)], size=len(genes))
        for g in genes:
            G.add_node(("G", g), kind="gene", comm=mg.loc[g, "Community"], bridge=bool(mg.loc[g, "Bridge"]))
            G.add_edge(("T", t), ("G", g))
    comms = sorted({d["comm"] for _, d in G.nodes(data=True)}, key=lambda c: int(c[1:]))
    col = {c: PAL[i % len(PAL)] for i, c in enumerate(comms)}
    pos = nx.spring_layout(G, seed=seed, k=0.9 / len(G) ** 0.5, iterations=300)
    nx.draw_networkx_edges(G, pos, ax=ax, width=0.25, alpha=0.25, edge_color="#777")
    genes = [n for n, d in G.nodes(data=True) if d["kind"] == "gene"]
    trms = [n for n, d in G.nodes(data=True) if d["kind"] == "term"]
    nx.draw_networkx_nodes(G, pos, nodelist=genes, ax=ax, node_size=14, node_shape="o",
                           node_color=[col[G.nodes[n]["comm"]] for n in genes],
                           edgecolors=["black" if G.nodes[n]["bridge"] else "none" for n in genes], linewidths=0.6)
    nx.draw_networkx_nodes(G, pos, nodelist=trms, ax=ax, node_shape="s",
                           node_size=[18 + 3 * G.nodes[n]["size"] for n in trms],
                           node_color=[col[G.nodes[n]["comm"]] for n in trms], edgecolors="white", linewidths=0.5, alpha=0.9)
    deg = mg["Degree"]
    top = list(deg.sort_values(ascending=False).index[:6])
    groups = []                                   # merge labels of genes drawn on top of each other
    for g in top:
        p = pos[("G", g)]
        for grp in groups:
            q = pos[("G", grp[0])]
            if ((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** 0.5 < 0.06:
                grp.append(g); break
        else:
            groups.append([g])
    for grp in groups:
        x, y = pos[("G", grp[0])]
        lab = (f"{grp[0]} ({deg[grp[0]]})" if len(grp) == 1 else
               ", ".join(grp) + f"\n(degree {min(deg[g] for g in grp)}–{max(deg[g] for g in grp)})")
        ax.annotate(lab, (x, y), xytext=(x + 0.12, y + 0.22) if len(grp) > 1 else (x + 0.02, y),
                    fontsize=6.5, ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#333") if len(grp) > 1 else None,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8) if len(grp) > 1 else None)
    ne = G.number_of_edges()
    ax.set_title(f"{title}\n{len(trms)} terms · {len(genes)} genes · {ne} edges · {len(comms)} communities · max gene degree {deg.max()}",
                 fontsize=9.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=col[c], label=f"{c} {NAMES[kind].get(c, '')}") for c in comms],
              loc="lower left", fontsize=7, frameon=False, ncol=1 if kind == "ora" else 4, handletextpad=0.1, columnspacing=0.8)
    ax.axis("off")

fig, axs = plt.subplots(1, 2, figsize=(14, 7))
draw(axs[0], "ora", "A. ORA, 10 collections (BH FDR < 0.05; each term with all its genes)")
draw(axs[1], "igea", "B. iGEA, 10 collections (q < 0.10; each gene at most once per collection)")
fig.text(0.5, 0.01, "Squares = terms (size = genes), circles = genes, black outline = bridge gene; colours = Leiden communities "
         "(term network, same procedure for both). Labels: six genes with the highest degree (number of terms).",
         ha="center", fontsize=7.5, color="#444")
plt.tight_layout(rect=(0, 0.03, 1, 1))
plt.savefig("out/ORA_vs_iGEA_10_collections_communities.png", dpi=200)
