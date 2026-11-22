"""Supplementary Figure 4: random gene lists vs the HIV network.
Usage: python supp_figure_04.py [<output folder>]   (written as sf5.py before the 21 Nov renumbering)"""
import sys as _sys, os as _os
_OUT = _sys.argv[1] if len(_sys.argv) > 1 else "."   # output folder
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, glob
sys.path.insert(0, str(_REPO) + "/code/ui"); sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from network_plot import _layout, _hull
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sf34lib import pull_in
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5})
PAL = ["#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD", "#17BECF", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22"]
split = lambda s: [g.strip() for g in str(s).split(",") if g.strip()]
R = str(_PUB) + "/analyses/"
RUNS = [(glob.glob(R + "step09_end_to_end/random_control_q1/random*/")[0], "A  Random list 1"),
        (glob.glob(R + "step12_supp_figures/random_seed1000/random*/")[0], "B  Random list 2"),
        (glob.glob(R + "step12_supp_figures/random_seed1001/random*/")[0], "C  Random list 3")]
def rep(d):
    l = open(glob.glob(d + "*statistical_report.txt")[0]).read().splitlines()
    n = l[1].split(": ")[1]; lc = l[2].split(": ", 1)[1]
    return n, lc
def draw(ax, d, title):
    t = pd.read_csv(d + "combined_iterative_results.tsv", sep="\t"); t = t[t["In network"] == True]
    terms = [dict(library=r.Library, term=r.Term, genes=split(r.Genes)) for r in t.itertuples()]
    mt = pd.read_csv(d + "network_modules_terms.tsv", sep="\t"); tc = {(r.Library, r.Term): r.Community for r in mt.itertuples()}
    mg = pd.read_csv(d + "network_modules_genes.tsv", sep="\t").set_index("Gene")
    pos = pull_in(_layout(terms, tc, {g: r for g, r in zip(mg.index, mg.itertuples())}, seed=2026), terms)
    comms = sorted(set(tc.values()), key=lambda c: int(c[1:])); col = {c: PAL[i % len(PAL)] for i, c in enumerate(comms)}
    for i, x in enumerate(terms):
        for g in set(x["genes"]):
            a, b = pos[("T", i)], pos[("G", g)]; ax.plot([a[0], b[0]], [a[1], b[1]], color="#999", lw=0.25, alpha=0.5, zorder=1)
    genes = sorted({g for x in terms for g in x["genes"]})
    ax.scatter([pos[("G", g)][0] for g in genes], [pos[("G", g)][1] for g in genes], s=4, zorder=2, c=[col[mg.loc[g, "Community"]] for g in genes], lw=0)
    ax.scatter([pos[("T", i)][0] for i in range(len(terms))], [pos[("T", i)][1] for i in range(len(terms))], marker="s", zorder=3,
               s=[8 + 1.2 * len(x["genes"]) for x in terms], c=[col[tc[(x["library"], x["term"])]] for x in terms], edgecolors="white", linewidths=0.3)
    n, lc = rep(d)
    epg = float(lc.split(", ")[4].split()[0])
    bb = ax.get_position(); ax.figure.text(bb.x0, 0.965, f"{title}\n{n} genes, {len(terms)} terms\nlargest cluster {lc.split(', ')[0]},\n{epg:.2f} edges per gene", fontsize=7.5, va="top")
    ax.set_aspect("equal"); ax.set_anchor("N"); ax.axis("off")
fig = plt.figure(figsize=(7.48, 5.2))
for k, (d, t) in enumerate(RUNS):
    draw(fig.add_axes([0.01 + k * 0.33, 0.44, 0.31, 0.42]), d, t)
null = pd.read_csv(glob.glob(R + "step09_end_to_end/HIV*/*random_lists.tsv")[0], sep="\t")
st = pd.read_csv(glob.glob(R + "step09_end_to_end/HIV*/*statistics.tsv")[0], sep="\t").set_index("Metric")
for k, (m, lab, key, bins) in enumerate([("cluster_size", "Nodes in the largest cluster", "Cluster Size", np.arange(0, 232, 4)),
                                          ("avg_edges_per_gene", "Edges per gene in the largest cluster", "Average Edges per Gene", np.arange(0, 4.6, 0.1))]):
    ax = fig.add_axes([0.08 + k * 0.48, 0.08, 0.40, 0.26])
    ax.hist(null[m], bins=bins, color="#9AA5B1", edgecolor="white", lw=0.3)
    obs = st.loc[key, "Observed"]; ax.axvline(obs, color="#D62728", lw=1.2)
    ax.set_yscale("log"); ax.set_xlabel(lab); ax.set_ylabel("Random lists")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(obs, ax.get_ylim()[1] * 0.6, f" HIV {obs:.3g}\n p = {st.loc[key, 'Empirical p']:.2g}" if st.loc[key, 'Empirical p'] > 0.0002 else f" HIV {obs:.3g}\n p < 0.0001",
            color="#D62728", fontsize=7, va="top", ha="right" if k == 0 else "left")
fig.text(0.01, 0.37, "D  HIV dependency factors (265 genes) vs 10,000 random lists of the same size, networks built after the term-level filter (q < 0.10)", fontsize=7.5)
fig.savefig(_OUT + "/Supp Figure 04 - Random lists vs HIV network.pdf"); fig.savefig(_OUT + "/Supp Figure 04 - Random lists vs HIV network.png", dpi=300)
