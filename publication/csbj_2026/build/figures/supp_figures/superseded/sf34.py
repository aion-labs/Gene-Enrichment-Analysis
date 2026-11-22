import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys
sys.path.insert(0, str(_REPO) + "/code/ui")
from network_plot import _layout, _hull
import pandas as pd, numpy as np, networkx as nx, matplotlib, textwrap
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family": "DejaVu Sans"})
PAL = ["#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD", "#17BECF", "#8C564B", "#E377C2"]
split = lambda s: [g.strip() for g in str(s).split(",") if g.strip()]
def nice(t):
    head, _, rest = t.partition(": ")
    return (head + ": " + rest.lower()) if rest else t.lower()
def pull_in(pos, terms):
    G = nx.Graph(); G.add_nodes_from(pos)
    G.add_edges_from((("T", i), ("G", g)) for i, x in enumerate(terms) for g in x["genes"])
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    P = {n: np.array(v, float) for n, v in pos.items()}
    main = np.array([P[n] for n in comps[0]]); c0 = main.mean(0); R = np.max(np.linalg.norm(main - c0, axis=1))
    for comp in comps[1:]:
        pts = np.array([P[n] for n in comp]); cc = pts.mean(0); r = np.max(np.linalg.norm(pts - cc, axis=1)) if len(pts) > 1 else 0
        d = cc - c0; d = d / (np.linalg.norm(d) or 1); new = c0 + d * (R + r + 0.12 * R)
        for n in comp: P[n] = P[n] - cc + new
    return {n: tuple(v) for n, v in P.items()}
def spread(pos, nt, it=400):
    keys = list(pos); P = np.array([pos[k] for k in keys], float)
    span = max(np.ptp(P[:, 0]), np.ptp(P[:, 1])); isT = np.array([k[0] == "T" for k in keys])
    for _ in range(it):
        D = P[:, None, :] - P[None, :, :]; dist = np.sqrt((D ** 2).sum(-1)) + 1e-9
        need = np.where(isT[:, None] & isT[None, :], 0.045, np.where(isT[:, None] | isT[None, :], 0.025, 0.012)) * span
        push = np.clip(need - dist, 0, None); np.fill_diagonal(push, 0)
        if push.max() < 1e-4 * span: break
        P += 0.5 * ((D / dist[..., None]) * push[..., None]).sum(1)
    return {k: tuple(P[j]) for j, k in enumerate(keys)}


def draw(d, out, size, title, gfs=5.5, tfs=6):
    from adjustText import adjust_text
    t = pd.read_csv(d + "combined_iterative_results.tsv", sep="\t"); t = t[t["In network"] == True]
    terms = [dict(library=r.Library, term=r.Term, genes=split(r.Genes)) for r in t.itertuples()]
    mt = pd.read_csv(d + "network_modules_terms.tsv", sep="\t"); tc = {(r.Library, r.Term): r.Community for r in mt.itertuples()}
    mg = pd.read_csv(d + "network_modules_genes.tsv", sep="\t").set_index("Gene")
    pos = pull_in(_layout(terms, tc, {g: r for g, r in zip(mg.index, mg.itertuples())}, seed=2026), terms)
    comms = sorted(set(tc.values()), key=lambda c: int(c[1:])); col = {c: PAL[i % len(PAL)] for i, c in enumerate(comms)}
    tcomm = {i: tc[(x["library"], x["term"])] for i, x in enumerate(terms)}
    pos = spread(pos, len(terms))
    # number terms by community, then by size
    order = sorted(range(len(terms)), key=lambda i: (int(tcomm[i][1:]), -len(terms[i]["genes"])))
    num = {i: k + 1 for k, i in enumerate(order)}
    fig = plt.figure(figsize=size); ax = fig.add_axes([0.0, 0.02, 0.66, 0.92])
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]; span = max(max(xs) - min(xs), max(ys) - min(ys))
    for c in comms:
        pts = [pos[("T", i)] for i, cc in tcomm.items() if cc == c] + [pos[("G", g)] for g in mg.index if mg.loc[g, "Community"] == c and not mg.loc[g, "Bridge"] and ("G", g) in pos]
        poly = _hull(pts, 0.035 * span); ax.fill([p[0] for p in poly], [p[1] for p in poly], color=col[c], alpha=0.08, lw=0)
    for i, x in enumerate(terms):
        for g in set(x["genes"]):
            a_, b_ = pos[("T", i)], pos[("G", g)]; ax.plot([a_[0], b_[0]], [a_[1], b_[1]], color="#999", lw=0.25, alpha=0.5, zorder=1)
    genes = sorted({g for x in terms for g in x["genes"]})
    ax.scatter([pos[("G", g)][0] for g in genes], [pos[("G", g)][1] for g in genes], s=5, zorder=2,
               c=[col[mg.loc[g, "Community"]] for g in genes], edgecolors=["black" if mg.loc[g, "Bridge"] else "none" for g in genes], linewidths=0.5)
    for i in range(len(terms)):
        p = pos[("T", i)]
        ax.text(p[0], p[1], str(num[i]), fontsize=tfs, ha="center", va="center", zorder=5, color="white", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc=col[tcomm[i]], ec="white", lw=0.4))
    texts = [ax.text(pos[("G", g)][0], pos[("G", g)][1], g, fontsize=gfs, zorder=4, color="#111",
                     fontweight="bold" if bool(mg.loc[g, "Bridge"]) else "normal") for g in genes]
    pad = 0.04 * span; ax.set_xlim(min(xs) - pad, max(xs) + pad); ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal"); ax.axis("off")
    adjust_text(texts, ax=ax, x=[pos[("G", g)][0] for g in genes] + [pos[("T", i)][0] for i in range(len(terms))],
                y=[pos[("G", g)][1] for g in genes] + [pos[("T", i)][1] for i in range(len(terms))],
                expand=(1.15, 1.3), arrowprops=dict(arrowstyle="-", color="#888", lw=0.3), max_move=None, time_lim=60)
    # key
    kx = fig.add_axes([0.67, 0.02, 0.32, 0.92]); kx.axis("off"); kx.set_xlim(0, 1)
    n_lines = len(terms) + 2 * len(comms); lh = 1.0 / max(n_lines, 60); y = 1.0
    mc = pd.read_csv(d + "network_modules_communities.tsv", sep="\t").set_index("Community")
    for c in comms:
        kx.text(0, y, f"{c}  ({mc.loc[c, 'Terms']} terms, {mc.loc[c, 'Genes']} genes)", fontsize=7, fontweight="bold", color=col[c], va="top"); y -= lh * 1.2
        for i in [i for i in order if tcomm[i] == c]:
            nm = nice(terms[i]["term"]); nm = nm if len(nm) <= 58 else nm[:56] + "…"
            kx.text(0.0, y, f"{num[i]:>2}", fontsize=6, va="top", ha="left", color=col[c], fontweight="bold", family="DejaVu Sans Mono")
            kx.text(0.07, y, f"{nm} ({len(terms[i]['genes'])})", fontsize=6, va="top"); y -= lh
        y -= lh * 0.6
    kx.set_ylim(min(y, 0), 1)
    fig.text(0.01, 0.985, title, fontsize=9, va="top")
    fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=300); print(out, len(terms), len(genes))
B = str(_PUB) + "/analyses/step09_end_to_end/"
G = str(_PUB) + "/analyses/step09_end_to_end/hiv_go_only/HIV_dependency_factors.symbols_20261120_114431/"
import os
if not os.path.exists(G + "network_modules_communities.tsv"):
    G = B + "hiv_go_only/HIV_dependency_factors.symbols_20261120_114431/"
draw(B + "HIV_dependency_factors.symbols_20261120_114347/", "figs/SF3", (11, 8.5),
     "HIV dependency factors, iGEA network with ten collections (56 terms, 168 genes, seven communities)")
draw(G, "figs/SF4", (11, 8.5), "HIV dependency factors, iGEA network with three Gene Ontology collections (Figure 2; 26 terms, 138 genes, six communities)", gfs=6, tfs=6.5)
