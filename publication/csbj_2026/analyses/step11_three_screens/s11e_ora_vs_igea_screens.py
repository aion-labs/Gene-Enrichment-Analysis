"""Same cross-screen comparison for ORA (FDR < 0.05, same gene sets/background/overlap as iGEA): term overlap and community
footprint matches, to compare with the iGEA result (s11b/s11c) on equal terms."""
import os, sys, glob, itertools
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "step10_extras"))
import common
from common import LIBS, BACKGROUND, FastCollection, genes_of
from network_modules import analyze_network_modules
H = os.path.dirname(os.path.abspath(__file__))
for s, f in [("Brass", "HIV_dependency_factors.symbols.txt"), ("König", "Konig_2008_HIV_dependency_factors.symbols.txt"), ("Zhou", "Zhou_2008_HIV_dependency_factors.symbols.txt")]:
    common.DATA[s] = f
cols = [FastCollection.from_library(l, BACKGROUND, 10, 600) for l in LIBS]
members = {(c.name, c.term_names[t]): {c.genes[j] for j in c.members[t]} for c in cols for t in range(len(c.term_names))}
allsets = list(members.values())
res = {}
for s in ["Brass", "König", "Zhou"]:
    g = set(genes_of(s)); terms = []
    for c in cols:
        x = np.zeros(len(c.genes)); [x.__setitem__(c.index[q], 1) for q in g if q in c.index]; k = c.A @ x
        p = np.asarray(c._pvalues("Fisher's Exact Test", k, c.K, x.sum()), dtype=float); q = multipletests(p, method="fdr_bh")[1]
        for t in np.where((q < 0.05) & (k >= 3))[0]:
            terms.append(dict(library=c.name, term=c.term_names[t], genes=[c.genes[j] for j in c.members[t] if c.genes[j] in g]))
    m = analyze_network_modules(terms)
    comm = {}
    for r in m.terms.itertuples(): comm.setdefault(r.Community, set()).add((r.Library, r.Term))
    deg = {}
    for t in terms:
        for q in t["genes"]: deg[q] = deg.get(q, 0) + 1
    res[s] = dict(terms={(t["library"], t["term"]) for t in terms}, comm=comm, n=len(terms), genes=len(deg),
                  edges=sum(deg.values()), mod=m.summary.get("modularity"), lead=dict(zip(m.communities.Community, m.communities["Lead term"])))
    print(s, "ORA terms", len(terms), "genes", len(deg), "edges/gene", round(sum(deg.values()) / max(1, len(deg)), 2), "communities", len(comm), "modularity", m.summary.get("modularity"))
rows = []
for a, b in itertools.combinations(res, 2):
    ta, tb = res[a]["terms"], res[b]["terms"]
    rows.append(dict(pair=f"{a}–{b}", ORA_terms_a=len(ta), ORA_terms_b=len(tb), shared=len(ta & tb), jaccard=round(len(ta & tb) / max(1, len(ta | tb)), 3)))
print(pd.DataFrame(rows).to_string(index=False))
print("ORA terms in all three:", len(res["Brass"]["terms"] & res["König"]["terms"] & res["Zhou"]["terms"]))
rng = np.random.default_rng(2026); cache = {}
def null(a, b):
    kk = tuple(sorted((a, b)))
    if kk not in cache:
        v = []
        for _ in range(2000):
            A = set().union(*[allsets[i] for i in rng.choice(len(allsets), a, replace=False)])
            B = set().union(*[allsets[i] for i in rng.choice(len(allsets), b, replace=False)])
            v.append(len(A & B) / len(A | B))
        cache[kk] = np.array(v)
    return cache[kk]
cm = []
for (s1, s2) in itertools.combinations(res, 2):
    for c1, t1 in res[s1]["comm"].items():
        for c2, t2 in res[s2]["comm"].items():
            A = set().union(*[members[t] for t in t1]); B = set().union(*[members[t] for t in t2]); j = len(A & B) / len(A | B)
            v = null(len(t1), len(t2)); p = (int((v >= j).sum()) + 1) / (len(v) + 1)
            cm.append(dict(screen_1=s1, community_1=c1, lead_1=res[s1]["lead"][c1], terms_1=len(t1), screen_2=s2, community_2=c2,
                           lead_2=res[s2]["lead"][c2], terms_2=len(t2), footprint_jaccard=round(j, 3), empirical_p=p))
cm = pd.DataFrame(cm).sort_values("footprint_jaccard", ascending=False); cm.to_csv(f"{H}/s11e_ora_community_similarity.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 45)
print(cm[cm.empirical_p <= 0.01].to_string(index=False))
print("pairs", len(cm), "p<=0.01:", int((cm.empirical_p <= 0.01).sum()))
for s in res: print(s, {c: (len(t), res[s]["lead"][c][:50]) for c, t in res[s]["comm"].items()})
