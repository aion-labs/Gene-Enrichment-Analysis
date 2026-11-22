import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5})
d = pd.read_csv(str(_PUB) + "/analyses/step12_supp_figures/s12a_size_trend.tsv", sep="\t").drop_duplicates()
d["q_filter"] = d["q_filter"].astype(str)
print(d.groupby(["q_filter", "list_size"]).size().unstack())
sizes = sorted(d.list_size.unique())
fig, axs = plt.subplots(2, 2, figsize=(7.48, 5.0), sharex=True)
for r, (q, lab) in enumerate([("none", "All selected terms (p < 0.01)"), ("0.10", "After the term-level filter (q < 0.10)")]):
    s = d[d.q_filter == q]
    for c, (m, yl) in enumerate([("cluster_size", "Nodes in the largest cluster"), ("avg_edges_per_gene", "Edges per gene")]):
        ax = axs[r, c]; data = [s[s.list_size == z][m].values for z in sizes]
        ax.boxplot(data, positions=range(len(sizes)), widths=0.6, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="#C9D3DE", lw=0.6), medianprops=dict(color="#1F4E79", lw=1.2), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        ax.plot(range(len(sizes)), [np.mean(x) for x in data], "o-", color="#D62728", ms=2.5, lw=0.8, label="mean")
        ax.set_ylabel(yl); ax.spines[["top", "right"]].set_visible(False)
        if c == 0: ax.set_title(("A  " if r == 0 else "B  ") + lab, loc="left", fontsize=8)
for ax in axs[1]: ax.set_xticks(range(len(sizes))); ax.set_xticklabels(sizes); ax.set_xlabel("Random gene list size")
axs[0, 1].legend(frameon=False, fontsize=7)
plt.tight_layout(); fig.savefig("figs/SF6.pdf"); fig.savefig("figs/SF6.png", dpi=300)
print(d.groupby(["q_filter","list_size"])[["cluster_size","avg_edges_per_gene"]].agg(["mean","median"]).round(2))
