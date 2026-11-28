"""Editor R1, point 3: multiple testing in the cross-screen community comparison.
Re-estimates the chance level of community footprint similarity with 50,000 random term groups per size pair
(s11c/s11e used 2,000, so no p-value could fall below 0.0005), for the 131 iGEA and the 16 ORA cross-screen
community pairs, and adjusts within each method (Benjamini-Hochberg FDR; Bonferroni shown for reference).
Random groups: terms drawn without replacement from all 10,439 analysed gene sets (10-600 genes after background
restriction), footprint = union of their genes; empirical p = (r + 1)/(n + 1).
Outputs: s11f_igea_pairs_adjusted.tsv, s11f_ora_pairs_adjusted.tsv, s11f_summary.tsv"""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "step10_extras"))
from common import LIBS, BACKGROUND, FastCollection
H = os.path.dirname(os.path.abspath(__file__))
N_DRAWS = int(os.environ.get("N_DRAWS", 50000)); SEED = 2026

sets = []
for l in LIBS:
    c = FastCollection.from_library(l, BACKGROUND, 10, 600)
    sets += [frozenset(c.genes[j] for j in c.members[t]) for t in range(len(c.term_names))]
genes = sorted(set().union(*sets)); gi = {g: i for i, g in enumerate(genes)}
W = (len(genes) + 63) // 64
M = np.zeros((len(sets), W), dtype=np.uint64)
for s, gs in enumerate(sets):
    for g in gs:
        i = gi[g]; M[s, i >> 6] |= np.uint64(1) << np.uint64(i & 63)
print(f"{len(sets)} gene sets, {len(genes)} genes, {N_DRAWS} draws per size pair")

rng = np.random.default_rng(SEED); cache = {}
def null(a, b):
    k = tuple(sorted((int(a), int(b))))
    if k not in cache:
        v = np.empty(N_DRAWS)
        for d in range(N_DRAWS):
            A = np.bitwise_or.reduce(M[rng.choice(len(sets), k[0], replace=False)], axis=0)
            B = np.bitwise_or.reduce(M[rng.choice(len(sets), k[1], replace=False)], axis=0)
            v[d] = np.bitwise_count(A & B).sum() / np.bitwise_count(A | B).sum()
        cache[k] = v
    return cache[k]

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

def run(df, size1, size2, exact_col="footprint_jaccard"):
    # same rule as s11c/s11e: random Jaccard >= observed Jaccard (rounded to 3 decimals)
    t0 = time.time(); ps, med = [], []
    for r in df.itertuples():
        v = null(size1(r), size2(r))
        ps.append((int((v >= getattr(r, exact_col)).sum()) + 1) / (N_DRAWS + 1)); med.append(round(float(np.median(v)), 3))
    df = df.copy(); df["random_median_jaccard"] = med; df["empirical_p"] = ps
    df["BH_FDR"] = bh(ps); df["bonferroni"] = np.minimum(np.array(ps) * len(ps), 1)
    print(f"  {len(df)} pairs in {time.time() - t0:.0f} s")
    return df.sort_values("empirical_p")

ig = pd.read_csv(f"{H}/s11b_community_similarity.tsv", sep="\t")
bm = pd.read_csv(f"{H}/s11b_community_best_matches.tsv", sep="\t")
size = {(r.screen, r.community): r.terms for r in bm.itertuples()}
ig = run(ig.drop(columns=["random_median_jaccard", "empirical_p"]),
         lambda r: size[(r.screen_1, r.community_1)], lambda r: size[(r.screen_2, r.community_2)])
ora = pd.read_csv(f"{H}/s11e_ora_community_similarity.tsv", sep="\t")
ora = run(ora.drop(columns=["empirical_p"]), lambda r: r.terms_1, lambda r: r.terms_2)
ig.to_csv(f"{H}/s11f_igea_pairs_adjusted.tsv", sep="\t", index=False)
ora.to_csv(f"{H}/s11f_ora_pairs_adjusted.tsv", sep="\t", index=False)
rows = []
for name, d in (("iGEA", ig), ("ORA", ora)):
    rows.append(dict(method=name, pairs=len(d), draws=N_DRAWS, p_le_0_01=int((d.empirical_p <= 0.01).sum()),
                     BH_lt_0_05=int((d.BH_FDR < 0.05).sum()), BH_lt_0_10=int((d.BH_FDR < 0.10).sum()),
                     bonferroni_lt_0_05=int((d.bonferroni < 0.05).sum())))
pd.DataFrame(rows).to_csv(f"{H}/s11f_summary.tsv", sep="\t", index=False)
np.savez_compressed(f"{H}/s11f_null_jaccard.npz", **{f"{a}_{b}": v for (a, b), v in cache.items()})  # random Jaccards per size pair
pd.set_option("display.width", 250)
print(pd.DataFrame(rows).to_string(index=False))
print(ig[ig.BH_FDR < 0.25].drop(columns=["shared_term_names"]).to_string(index=False))
print(ora[ora.BH_FDR < 0.25].to_string(index=False))
