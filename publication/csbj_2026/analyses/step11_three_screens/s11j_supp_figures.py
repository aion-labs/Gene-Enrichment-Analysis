"""Supplementary figures for the three-screen comparison (R2, editor point 3).
(1) s11j_brass_ranked_bars: each Brass community vs all König and Zhou communities, ranked by footprint Jaccard, with the
    median Jaccard of matched random groups (chance level) and the adjusted significance (matched null, K = 25, BH over 131 pairs).
(2) s11j_multiple_testing: (a) observed vs matched-null median Jaccard for all 131 iGEA and 16 ORA pairs; (b) distribution of the
    131 iGEA p-values; (c) FDR of the leading pairs for K = 10, 25, 50 nearest gene sets in size.
Inputs: s11i_footprint_jaccard_matrix.tsv, s11h_igea_pairs_matched{,_K10,_K50}.tsv, s11h_ora_pairs_matched{...}.tsv"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
H = os.path.dirname(os.path.abspath(__file__))
COL = {"Brass": "#1f77b4", "König": "#ff7f0e", "Zhou": "#2ca02c"}
M = pd.read_csv(os.path.join(H, "s11i_footprint_jaccard_matrix.tsv"), sep="\t", index_col=0)
lab = list(M.index); key = {l: tuple(l.split(" · ")[0].split(":")) for l in lab}
T = {k: pd.read_csv(os.path.join(H, f"s11h_{m}_pairs_matched{k}.tsv"), sep="\t") for m in ("igea",) for k in ("", "_K10", "_K50")}
O = {k: pd.read_csv(os.path.join(H, f"s11h_ora_pairs_matched{k}.tsv"), sep="\t") for k in ("", "_K10", "_K50")}
def lut(d, col):
    out = {}
    for r in d.itertuples():
        out[(r.screen_1, r.community_1, r.screen_2, r.community_2)] = out[(r.screen_2, r.community_2, r.screen_1, r.community_1)] = getattr(r, col)
    return out
FDR, MED = lut(T[""], "matched_BH_FDR"), lut(T[""], "matched_median_jaccard")
def star(q): return "*" if q < 0.05 else "†" if q < 0.10 else ""

# (1) ranked bars
q = [l for l in lab if l.startswith("Brass:")]
fig, axes = plt.subplots(4, 2, figsize=(15, 19)); axes = axes.ravel()
for ax, a in zip(axes, q):
    others = [l for l in lab if not l.startswith("Brass:")]
    vals = sorted(((o, M.loc[a, o]) for o in others), key=lambda x: -x[1])
    y = np.arange(len(vals))[::-1]
    for yi, (o, v) in zip(y, vals):
        s = key[o][0]; k = (*key[a], *key[o]); f = FDR[k]
        ax.barh(yi, v, color=COL[s] if s == "König" else "white", edgecolor=COL[s], hatch=None if s == "König" else "////", lw=1.4, height=0.72)
        ax.plot([MED[k]] * 2, [yi - 0.38, yi + 0.38], color="black", lw=1.6)
        ax.text(max(v, MED[k]) + 0.012, yi, f"{v:.2f}{star(f)}", va="center", fontsize=8.5, fontweight="bold" if f < 0.10 else "normal")
        ax.text(-0.012, yi, o, ha="right", va="center", fontsize=8.3, color=COL[s], fontweight="bold" if f < 0.10 else "normal")
    ax.set_yticks([]); ax.set_xlim(0, 0.9); ax.set_ylim(-0.6, len(vals) - 0.4)
    ax.set_title(a, loc="left", fontsize=11, fontweight="bold", color=COL["Brass"])
    ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(labelsize=8); ax.set_xlabel("footprint Jaccard", fontsize=8.5)
ax = axes[-1]; ax.axis("off")
ax.legend(handles=[Patch(facecolor=COL["König"], edgecolor=COL["König"], label="König community"),
                   Patch(facecolor="white", edgecolor=COL["Zhou"], hatch="////", label="Zhou community"),
                   Line2D([0], [0], color="black", lw=1.6, label="chance level: median Jaccard of\nmatched random groups of gene sets")],
          loc="center", frameon=False, fontsize=11, title="* FDR < 0.05, † FDR < 0.10 (bold)\nmatched null, BH over 131 pairs", title_fontsize=11)
fig.suptitle("Each community of the Brass network compared with every community of the König and Zhou networks",
             fontsize=13, fontweight="bold", y=0.995)
fig.tight_layout(rect=(0.2, 0, 1, 0.98)); fig.subplots_adjust(wspace=1.05)
out = os.path.join(H, "s11j_brass_ranked_bars"); fig.savefig(out + ".pdf", bbox_inches="tight"); fig.savefig(out + ".png", dpi=200, bbox_inches="tight")

# (2) multiple testing overview
fig = plt.figure(figsize=(19, 6.6)); gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 0.75, 0.62, 1.0], wspace=0.3)
ax = fig.add_subplot(gs[0]); d = T[""]
cat = np.where(d.matched_BH_FDR < 0.05, 0, np.where(d.matched_BH_FDR < 0.10, 1, 2))
for c, (lbl, colr) in enumerate([("iGEA, FDR < 0.05", "#d62728"), ("iGEA, FDR < 0.10", "#ff9f40"), ("iGEA, not significant", "#b0b0b0")]):
    x = d[cat == c]; ax.scatter(x.matched_median_jaccard, x.footprint_jaccard, s=34 if c < 2 else 18, c=colr, edgecolor="black" if c < 2 else "none", lw=0.5, label=f"{lbl} ({len(x)})", zorder=3)
o = O[""]; ax.scatter(o.matched_median_jaccard, o.footprint_jaccard, marker="^", s=46, facecolor="white", edgecolor="black", lw=1, label=f"ORA, not significant ({len(o)})", zorder=4)
ax.plot([0, 0.8], [0, 0.8], ls="--", color="#666", lw=1, label="observed = chance")
ax.set_xlim(0, 0.25); ax.set_ylim(0, 0.82); ax.set_xlabel("chance level: median Jaccard of matched random groups"); ax.set_ylabel("observed footprint Jaccard")
ax.legend(frameon=False, fontsize=8.5, loc="upper right"); ax.set_title("a  Observed overlap vs chance, all community pairs", loc="left", fontweight="bold", fontsize=11)
for r in d[d.matched_BH_FDR < 0.10].itertuples():
    t = f"{r.screen_1[0]}{r.community_1[1:]}–{r.screen_2[0]}{r.community_2[1:]}"
    off = {"B7–Z5": (-12, 7), "K2–Z5": (5, -9), "K5–Z2": (6, -3), "B5–Z2": (6, 3)}.get(t, (5, 2))
    ax.annotate(t, (r.matched_median_jaccard, r.footprint_jaccard), fontsize=8, xytext=off, textcoords="offset points")
ax = fig.add_subplot(gs[1])
ax.hist(d.matched_p, bins=np.linspace(0, 1, 21), color="#9ab", edgecolor="white")
ax.axhline(len(d) / 20, color="black", ls="--", lw=1); ax.text(0.98, len(d) / 20 + 0.6, "uniform", ha="right", fontsize=8)
ax.set_xlabel("matched-null p-value"); ax.set_ylabel("iGEA community pairs (of 131)")
ax.set_title("b  p-values of the 131 iGEA pairs", loc="left", fontweight="bold", fontsize=11)
ax = fig.add_subplot(gs[3])
F = {k: lut(T[k], "matched_BH_FDR") for k in T}
lead = d[(d.matched_BH_FDR < 0.25) | (T["_K10"].set_index(["screen_1", "community_1", "screen_2", "community_2"]).loc[list(zip(d.screen_1, d.community_1, d.screen_2, d.community_2)), "matched_BH_FDR"].values < 0.10)]
SH = {("Brass", "C3"): "Mediator", ("König", "C3"): "Mediator", ("Zhou", "C3"): "Mediator", ("Brass", "C7"): "cullin-RING", ("König", "C2"): "proteasome",
      ("Zhou", "C5"): "K48 ubiq.", ("Brass", "C4"): "nuclear pore", ("König", "C4"): "nuclear pore", ("König", "C5"): "RNA pol II", ("Zhou", "C2"): "RNA pol II",
      ("Brass", "C5"): "RNA processing", ("König", "C1"): "mRNA splicing", ("Brass", "C1"): "Golgi", ("Brass", "C2"): "AKT/NF-κB", ("Zhou", "C1"): "NF-κB/AKT"}
rows = [(f"{r.screen_1}:{r.community_1} {SH.get((r.screen_1, r.community_1), '')} – {r.screen_2}:{r.community_2} {SH.get((r.screen_2, r.community_2), '')}",
         [F[k][(r.screen_1, r.community_1, r.screen_2, r.community_2)] for k in ("_K10", "", "_K50")]) for r in lead.itertuples()]
rows.sort(key=lambda x: x[1][1])
for i, (nm, v) in enumerate(rows[::-1]):
    for x, mk, c in zip(v, ("s", "o", "D"), ("#6a3d9a", "black", "#1b9e77")):
        ax.scatter(max(x, 1e-4), i, marker=mk, s=30 if mk == "o" else 22, color=c, zorder=3)
    ax.plot([max(min(v), 1e-4), max(v)], [i, i], color="#bbb", lw=1, zorder=1)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=8)
ax.set_xscale("log"); ax.axvline(0.05, color="#d62728", ls="--", lw=1); ax.axvline(0.10, color="#ff9f40", ls=":", lw=1)
ax.set_xlabel("Benjamini–Hochberg FDR (131 pairs)")
ax.legend(handles=[Line2D([0], [0], marker=m, ls="", color=c, label=l) for m, c, l in (("s", "#6a3d9a", "K = 10"), ("o", "black", "K = 25 (main)"), ("D", "#1b9e77", "K = 50"))],
          frameon=True, framealpha=0.9, fontsize=8.5, loc="upper right")
ax.set_title("c  Sensitivity to the size matching (K)", loc="left", fontweight="bold", fontsize=11)
fig.suptitle("Multiple testing of the cross-screen community comparison (annotation-matched null, 50,000 random pairs per community pair)",
             fontsize=12.5, fontweight="bold", y=1.02)
out = os.path.join(H, "s11j_multiple_testing"); fig.savefig(out + ".pdf", bbox_inches="tight"); fig.savefig(out + ".png", dpi=200, bbox_inches="tight")
print("ok", len(rows), "pairs in panel c")
