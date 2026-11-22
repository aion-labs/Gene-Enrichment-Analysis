"""Supplementary Figure 1: ORA vs iGEA terms per collection (from Supplementary Tables 4 and 5).
Usage: python supp_figure_01.py [<output folder>]   (written as sf1.py before the 21 Nov renumbering)"""
import sys as _sys, os as _os
_OUT = _sys.argv[1] if len(_sys.argv) > 1 else "."   # output folder
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
D = str(_PUB) + "/supplement/"
ora = pd.read_csv(D + "Supp Table 04 - ORA results.txt", sep="\t"); ig = pd.read_csv(D + "Supp Table 05 - iGEA results.txt", sep="\t")
ig["kept"] = ig["In network"].astype(str) == "True"
order = ["H: Hallmark Gene Sets", "C2: BioCarta", "C2: KEGG Legacy", "C2: KEGG MEDICUS", "C2: Pathway Interaction Database",
         "C2: Reactome Pathways", "C2: WikiPathways", "C5: Gene Ontology: Biological Process", "C5: Gene Ontology: Cellular Component",
         "C5: Gene Ontology: Molecular Function"]
assert set(ora.Collection) <= set(order) and set(ig.Collection) <= set(order), (set(ora.Collection) | set(ig.Collection)) - set(order)
short = {o: o.split(": ", 1)[1].replace("Gene Ontology: ", "GO ").replace(" Gene Sets", "").replace(" Pathways", "") for o in order}
def counts(o_sel, i_sel):
    rows = []
    for c in order:
        a = set(ora[(ora.Collection == c) & o_sel].Term); b = set(ig[(ig.Collection == c) & i_sel].Term)
        rows.append((len(a & b), len(b - a), len(a - b)))
    return rows
P = [("A  ORA p < 0.01 and iGEA p < 0.01", counts(ora["p-value"] < 0.01, ig["Iteration p-value"] < 0.01)),
     ("B  ORA FDR < 0.05 and iGEA q < 0.10", counts(ora["FDR"] < 0.05, ig["kept"]))]
C = ["#3B6EA8", "#2BA59A", "#E08A2E"]; L = ["Both", "iGEA only", "ORA only"]
fig, axs = plt.subplots(1, 2, figsize=(7.48, 3.5), sharey=True)
for ax, (title, rows) in zip(axs, P):
    y = range(len(order))[::-1]; left = [0] * len(order)
    for k in range(3):
        v = [r[k] for r in rows]
        ax.barh(list(y), v, left=left, color=C[k], label=L[k], height=0.7)
        left = [a + b for a, b in zip(left, v)]
    for yy, tot, r in zip(y, left, rows):
        ax.text(tot + 1.5, yy, f"{r[0]}/{r[1]}/{r[2]}", va="center", fontsize=7, color="#444")
    tot = [sum(x) for x in zip(*rows)]
    ax.set_title(f"{title}\nboth {tot[0]}, iGEA only {tot[1]}, ORA only {tot[2]}", fontsize=8, loc="left"); ax.set_xlabel("Terms")
    ax.spines[["top", "right"]].set_visible(False); ax.set_xlim(0, max(left) * 1.25)
axs[0].set_yticks(list(range(len(order))[::-1])); axs[0].set_yticklabels([short[o] for o in order])
fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, fontsize=7.5, loc="lower center", ncol=3)
plt.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(_OUT + "/Supp Figure 01 - ORA vs iGEA terms per collection.pdf"); fig.savefig(_OUT + "/Supp Figure 01 - ORA vs iGEA terms per collection.png", dpi=300)
for t, rows in P: print(t, [sum(x) for x in zip(*rows)], rows)
