"""Worked example of one cross-screen community comparison (editor R1, point 3): Brass C4 vs König C4 (nuclear pore).
(a) Jaccard index between the full gene sets of the terms of the two communities.
(b) Screen hits: input genes of each screen that belong to the terms of the two communities.
(c) Footprints (Venn): union of the full gene sets of each community's terms; Jaccard of the two footprints.
(d) Chance level: footprint Jaccard of the 50,000 pairs of random term groups of the same sizes used in s11f.
Output: s11g_example_pair_nuclear_pore.pdf / .png"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
H = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(H, "s11b_compare.py")).read().split("rows = []")[0]
__file__ = os.path.join(H, "s11b_compare.py"); exec(src)       # loads members, data (communities) and inputs (screen hits)
S1, C1, S2, C2 = "Brass", "C4", "König", "C4"
T1, T2 = sorted(data[S1]["comm"][C1], key=lambda t: -len(members[t])), sorted(data[S2]["comm"][C2], key=lambda t: -len(members[t]))
short = lambda t: t[1].split("_", 1)[1].replace("_", " ").lower().replace("rna", "RNA").replace("mrna", "mRNA").replace("gtpase", "GTPase").replace("npc", "NPC") + f"  [{t[1].split('_')[0]}]"
J = np.array([[len(members[a] & members[b]) / len(members[a] | members[b]) for b in T2] for a in T1])
F1 = set().union(*[members[t] for t in T1]); F2 = set().union(*[members[t] for t in T2])
jf = len(F1 & F2) / len(F1 | F2)
hits1, hits2 = inputs[S1], inputs[S2]
adj = pd.read_csv(os.path.join(H, "s11f_igea_pairs_adjusted.tsv"), sep="\t")
row = adj[(adj.screen_1 == S1) & (adj.community_1 == C1) & (adj.screen_2 == S2) & (adj.community_2 == C2)].iloc[0]

# chance level: the random Jaccards of the main analysis (s11f, 50,000 draws for this pair of community sizes)
nz = np.load(os.path.join(H, "s11f_null_jaccard.npz"))
null = nz["_".join(map(str, sorted((len(T1), len(T2)))))]; n = len(null)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
CB, CK, CBK = "#1f77b4", "#ff7f0e", "#7f3fbf"          # Brass, König, both
fig = plt.figure(figsize=(15, 18))
gs = fig.add_gridspec(3, 2, width_ratios=[1.25, 1], height_ratios=[1, 1, 1.45], wspace=0.35, hspace=1.05)

def screen_tag(ax, text, color, x, y, rot=0, ha="center", va="center"):
    ax.text(x, y, text, transform=ax.transAxes, rotation=rot, ha=ha, va=va, fontsize=10, fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.35", fc=color, ec="none"), clip_on=False)

# (a) term x term Jaccard heatmap
ax = fig.add_subplot(gs[0:2, 0])
im = ax.imshow(J, cmap="Greys", vmin=0, vmax=1, aspect="auto")
for i in range(len(T1)):
    for j in range(len(T2)):
        ax.text(j, i, f"{J[i, j]:.2f}", ha="center", va="center", fontsize=7, color="white" if J[i, j] > 0.55 else "black")
        if T1[i] == T2[j]:
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#d62728", lw=2))
ax.set_yticks(range(len(T1)), [short(t) for t in T1]); ax.set_xticks(range(len(T2)), [short(t) for t in T2], rotation=40, ha="right")
for lab in ax.get_yticklabels(): lab.set_color(CB)
for lab in ax.get_xticklabels(): lab.set_color(CK)
# the two screens on the two axes: thick colored axis lines and a named tag on each axis
for side, c in (("left", CB), ("bottom", CK)):
    ax.spines[side].set_color(c); ax.spines[side].set_linewidth(5)
ax.tick_params(axis="y", colors=CB, length=6, width=2); ax.tick_params(axis="x", colors=CK, length=6, width=2)
screen_tag(ax, f"ROWS: BRASS screen · community {C1} ({len(T1)} terms)", CB, -0.95, 0.5, rot=90)
screen_tag(ax, f"COLUMNS: KÖNIG screen · community {C2} ({len(T2)} terms)", CK, 0.5, 1.035, va="bottom")
ax.set_title("a  Jaccard index between the gene sets of the terms\n(red box: the same term was kept in both screens)",
             loc="left", fontsize=10, fontweight="bold", pad=46)
fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="Jaccard (gene sets)")

# (c) footprints as a Venn diagram
from matplotlib.patches import Circle
ax = fig.add_subplot(gs[0, 1]); ax.set_aspect("equal"); ax.axis("off")
o1, o2, sh = len(F1 - F2), len(F2 - F1), len(F1 & F2)
ax.add_patch(Circle((0.36, 0.5), 0.3, fc=CB, alpha=0.35, ec=CB, lw=2.5))
ax.add_patch(Circle((0.64, 0.5), 0.3, fc=CK, alpha=0.35, ec=CK, lw=2.5))
ax.text(0.2, 0.5, f"{o1}", ha="center", va="center", fontsize=13, fontweight="bold")
ax.text(0.5, 0.5, f"{sh}", ha="center", va="center", fontsize=13, fontweight="bold")
ax.text(0.8, 0.5, f"{o2}", ha="center", va="center", fontsize=13, fontweight="bold")
ax.text(0.2, 0.40, "Brass only", ha="center", fontsize=8); ax.text(0.5, 0.40, "both", ha="center", fontsize=8)
ax.text(0.8, 0.40, "König only", ha="center", fontsize=8)
screen_tag(ax, f"BRASS {C1}\nfootprint\n{len(F1)} genes", CB, -0.02, 0.5, ha="right")
screen_tag(ax, f"KÖNIG {C2}\nfootprint\n{len(F2)} genes", CK, 1.02, 0.5, ha="left")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.05)
ax.set_title(f"c  Community footprints (union of the full gene sets of each community's terms)\n"
             f"Jaccard = shared / union = {sh} / {len(F1 | F2)} = {jf:.3f}", loc="left", fontsize=10, fontweight="bold")
h1 = len(F1 & F2 & hits1); h2 = len(F1 & F2 & hits2); hb = len(F1 & F2 & hits1 & hits2)
ax.text(0.0, 0.08, f"Screen hits inside the {sh} shared genes: Brass {h1}, König {h2}, hits in both screens {hb}.\n"
        "The footprints overlap because both communities describe the nuclear pore and RNA export,\n"
        "not because the screens found the same genes.", fontsize=8, va="top")

# (b) screen hits x terms: groups labeled in text and cells patterned, so color is not needed to read it
ax = fig.add_subplot(gs[2, :])
terms = [t for t in T1 if t not in T2] + [t for t in T1 if t in T2] + [t for t in T2 if t not in T1]
g1 = sorted((hits1 - hits2) & set().union(*[members[t] for t in T1]))
g2 = sorted((hits2 - hits1) & set().union(*[members[t] for t in T2]))
gb = sorted(hits1 & hits2 & (set().union(*[members[t] for t in T1]) | set().union(*[members[t] for t in T2])))
genes_b = g1 + gb + g2
style = lambda g: (CBK, "xx") if g in gb else ((CB, None) if g in g1 else (CK, "////"))
for i, g in enumerate(genes_b):
    c, h = style(g)
    for j, t in enumerate(terms):
        if g in members[t]:
            ax.add_patch(Rectangle((i - 0.45, j - 0.45), 0.9, 0.9, fc=c, ec="white" if h else c, hatch=h, lw=0))
        else:
            ax.add_patch(Rectangle((i - 0.45, j - 0.45), 0.9, 0.9, color="#f2f2f2"))
ax.set_xlim(-0.6, len(genes_b) - 0.4); ax.set_ylim(len(terms) - 0.4, -0.6)
ax.set_xticks(range(len(genes_b)), genes_b, rotation=90, fontsize=6.5)
for lab, g in zip(ax.get_xticklabels(), genes_b):
    lab.set_color(style(g)[0])
n1, nsh = len([t for t in T1 if t not in T2]), len([t for t in T1 if t in T2])
ax.set_yticks(range(len(terms)), [short(t) for t in terms])
for lab, t in zip(ax.get_yticklabels(), terms):
    lab.set_color(CBK if (t in T1 and t in T2) else (CB if t in T1 else CK))
for y in (n1 - 0.5, n1 + nsh - 0.5):
    ax.axhline(y, color="black", lw=1.2)
for x in (len(g1) - 0.5, len(g1) + len(gb) - 0.5):
    ax.axvline(x, color="black", lw=1.2)
# column group headers (genes) and row group labels (terms)
for x0, x1, txt, c in ((0, len(g1), f"BRASS hits only ({len(g1)} genes)\nsolid", CB),
                       (len(g1), len(g1) + len(gb), f"hits in\nboth ({len(gb)})\ncross-hatched", CBK),
                       (len(g1) + len(gb), len(genes_b), f"KÖNIG hits only ({len(g2)} genes)\nhatched", CK)):
    ax.annotate("", xy=(x0 - 0.45, -1.0), xytext=(x1 - 0.55, -1.0), xycoords="data",
                arrowprops=dict(arrowstyle="-", color=c, lw=3), annotation_clip=False)
    ax.text((x0 + x1 - 1) / 2, -1.35, txt, ha="center", va="bottom", fontsize=9, fontweight="bold", color=c)
for y0, y1, txt, c in ((0, n1, "terms of the\nBRASS community\nonly", CB), (n1, n1 + nsh, "terms kept in\nBOTH communities", CBK),
                       (n1 + nsh, len(terms), "terms of the\nKÖNIG community\nonly", CK)):
    ax.annotate("", xy=(len(genes_b) + 0.1, y0 - 0.45), xytext=(len(genes_b) + 0.1, y1 - 0.55), xycoords="data",
                arrowprops=dict(arrowstyle="-", color=c, lw=3), annotation_clip=False)
    ax.text(len(genes_b) + 0.5, (y0 + y1 - 1) / 2, txt, ha="left", va="center", fontsize=9, fontweight="bold", color=c)
ax.set_title("b  Screen hits (columns) in the gene sets of the terms of both communities (rows); filled = the gene belongs to the term's gene set\n\n\n\n",
             loc="left", fontsize=10, fontweight="bold")
for s_ in ax.spines.values(): s_.set_visible(False)

# (d) chance level
ax = fig.add_subplot(gs[1, 1])
ax.hist(null, bins=np.linspace(0, max(0.4, null.max()), 81), color="#9aa5b1")
ax.axvline(jf, color="#d62728", lw=2)
r = int((null >= round(jf, 3)).sum())
ax.text(jf, ax.get_ylim()[1] * 0.9, f"  observed {jf:.3f}", color="#d62728", fontsize=9, va="top")
ax.set_yscale("log"); ax.set_xlabel("footprint Jaccard of two random term groups"); ax.set_ylabel("random pairs (log scale)")
ax.set_title(f"d  Chance level: {n:,} pairs of random groups\nof {len(T1)} and {len(T2)} terms (median Jaccard {np.median(null):.3f})",
             loc="left", fontsize=10, fontweight="bold")
ax.text(0.98, 0.55, f"empirical p = (r + 1)/(n + 1) = {(r + 1) / (n + 1):.5f}\n"
        f"{r} of {n:,} random pairs reach {jf:.3f}\n"
        f"Benjamini–Hochberg over the 131 iGEA pairs: FDR = {row.BH_FDR:.3f}", transform=ax.transAxes, ha="right", fontsize=8.5,
        bbox=dict(fc="white", ec="#cccccc"))
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
fig.suptitle(f"How two communities from different screens are compared: {S1} {C1} and {S2} {C2} (nuclear pore and RNA export)",
             fontsize=12, fontweight="bold", x=0.02, ha="left")
out = os.path.join(H, "s11g_example_pair_nuclear_pore")
fig.savefig(out + ".pdf", bbox_inches="tight"); fig.savefig(out + ".png", dpi=200, bbox_inches="tight")
print("saved", out, "| J", round(jf, 3), "| r", r, "| hits", len(g1), len(g2), len(gb))
