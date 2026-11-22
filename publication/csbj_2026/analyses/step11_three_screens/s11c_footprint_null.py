"""Chance level of community footprint similarity: for each cross-screen community pair, Jaccard of the footprints of two
random term groups of the same sizes (terms drawn from all 10,439 analysed gene sets, 2,000 draws). Empirical p = (r+1)/(n+1)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "step10_extras"))
from common import LIBS, BACKGROUND, FastCollection
H = os.path.dirname(os.path.abspath(__file__))
sets = []
for l in LIBS:
    c = FastCollection.from_library(l, BACKGROUND, 10, 600)
    sets += [frozenset(c.genes[j] for j in c.members[t]) for t in range(len(c.term_names))]
cm = pd.read_csv(f"{H}/s11b_community_similarity.tsv", sep="\t")
bm = pd.read_csv(f"{H}/s11b_community_best_matches.tsv", sep="\t")
size = {(r.screen, r.community): r.terms for r in bm.itertuples()}
rng = np.random.default_rng(2026); cache = {}
def null(a, b):
    k = tuple(sorted((a, b)))
    if k not in cache:
        v = []
        for _ in range(2000):
            A = set().union(*[sets[i] for i in rng.choice(len(sets), a, replace=False)])
            B = set().union(*[sets[i] for i in rng.choice(len(sets), b, replace=False)])
            v.append(len(A & B) / len(A | B))
        cache[k] = np.array(v)
    return cache[k]
ps, med = [], []
for r in cm.itertuples():
    v = null(size[(r.screen_1, r.community_1)], size[(r.screen_2, r.community_2)])
    ps.append((int((v >= r.footprint_jaccard).sum()) + 1) / (len(v) + 1)); med.append(round(float(np.median(v)), 3))
cm["random_median_jaccard"] = med; cm["empirical_p"] = ps
cm.to_csv(f"{H}/s11b_community_similarity.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(cm[cm.empirical_p < 0.05].drop(columns=["shared_term_names"]).to_string(index=False))
print("pairs tested:", len(cm), "p<0.05:", (cm.empirical_p < 0.05).sum(), "p<0.001:", (cm.empirical_p < 0.001).sum())
