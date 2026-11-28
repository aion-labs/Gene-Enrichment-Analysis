"""Editor R1, point 3 (refined null): chance level of community footprint similarity with random term groups MATCHED to
each community, term by term. Each term of a community is replaced by a gene set drawn from
  - the same collection,
  - among all analysed gene sets of that collection (10-600 genes after background restriction), with no reference to
    any screen's input genes,
  - of similar size: the K gene sets of that collection closest in size to the term (default K = 25; the term itself excluded).
Within a random group no gene set is used twice.
Footprint = union of the full gene sets of a group; similarity = Jaccard of the two footprints; empirical
p = (r + 1)/(n + 1) over N = 50,000 independent pairs of matched random groups; Benjamini-Hochberg within each method
(131 iGEA pairs, 16 ORA pairs). The uniform null of s11f (terms drawn from all gene sets) is kept for comparison.
Outputs (TAG = '' for K = 25, '_K10' / '_K50' for the sensitivity runs): s11h_igea_pairs_matched{TAG}.tsv,
s11h_ora_pairs_matched{TAG}.tsv, s11h_summary{TAG}.tsv, s11h_example_null{TAG}.npz"""
import os, sys, time, itertools
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
N, K, CHUNK, SEED = int(os.environ.get("N_DRAWS", 50000)), int(os.environ.get("K", 25)), 2000, 2026
TAG = os.environ.get("TAG", "")

def load(script, stop):
    ns = {"__file__": os.path.join(H, script), "__name__": "x"}
    exec(open(os.path.join(H, script)).read().split(stop)[0], ns); return ns
ig = load("s11b_compare.py", "rows = []")          # members, data[screen]["comm"], inputs, LIBS, BACKGROUND, FastCollection
ora = load("s11e_ora_vs_igea_screens.py", "rows = []")   # res[screen]["comm"] (ORA communities, FDR < 0.05)
members, inputs = ig["members"], ig["inputs"]
for k in ora["members"]:
    assert k in members, k                      # same gene sets and keys in both scripts
keys = list(members); kidx = {k: i for i, k in enumerate(keys)}
genes = sorted(set().union(*members.values())); gi = {g: i for i, g in enumerate(genes)}; W = (len(genes) + 63) // 64
M = np.zeros((len(keys), W), dtype=np.uint64)
for s, k in enumerate(keys):
    for g in members[k]:
        i = gi[g]; M[s, i >> 6] |= np.uint64(1) << np.uint64(i & 63)
size = np.array([len(members[k]) for k in keys]); coll = np.array([k[0] for k in keys])
print(f"{len(keys)} gene sets; K = {K} size-matched candidates per term", flush=True)

def candidates(screen, term):
    pool = np.where((coll == term[0]) & (np.arange(len(keys)) != kidx[term]))[0]
    d = np.abs(size[pool] - len(members[term]))
    return pool[np.argsort(d, kind="stable")[:K]]

def footprints(cands, n, rng):
    """n random groups (one draw per term from its candidates) -> (n, W) footprint bitsets"""
    out = np.empty((n, W), dtype=np.uint64)
    for a in range(0, n, CHUNK):
        b = min(n, a + CHUNK)
        idx = np.stack([rng.choice(c, b - a) for c in cands], axis=1)
        for _ in range(1000):                                   # no gene set twice within a group: redraw clashes
            srt = np.sort(idx, axis=1); dup = np.where((srt[:, 1:] == srt[:, :-1]).any(axis=1))[0]
            if not len(dup): break
            for r in dup:
                seen = set()
                for j, c in enumerate(cands):
                    if idx[r, j] in seen:
                        free = [x for x in c if x not in seen]
                        idx[r, j] = rng.choice(free) if free else idx[r, j]
                    seen.add(idx[r, j])
        out[a:b] = np.bitwise_or.reduce(M[idx], axis=1)
    return out

def jac(A, B):
    return np.bitwise_count(A & B).sum(1) / np.bitwise_count(A | B).sum(1)

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

def run(comms, label, example=None):
    rng = np.random.default_rng(SEED); rows = []; t0 = time.time(); ex = None
    names = [(s, c) for s in comms for c in sorted(comms[s], key=lambda x: int(x[1:]))]
    cands = {(s, c): [candidates(s, t) for t in sorted(comms[s][c])] for s, c in names}
    fp = {(s, c): set().union(*[members[t] for t in comms[s][c]]) for s, c in names}
    for (s1, c1), (s2, c2) in itertools.combinations(names, 2):
        if s1 == s2: continue
        A, B = fp[(s1, c1)], fp[(s2, c2)]; j = round(len(A & B) / len(A | B), 3)
        v = jac(footprints(cands[(s1, c1)], N, rng), footprints(cands[(s2, c2)], N, rng))
        r = int((v >= j).sum())
        rows.append(dict(screen_1=s1, community_1=c1, terms_1=len(comms[s1][c1]), screen_2=s2, community_2=c2,
                         terms_2=len(comms[s2][c2]), footprint_jaccard=j, matched_median_jaccard=round(float(np.median(v)), 3),
                         matched_p=(r + 1) / (N + 1)))
        if example == ((s1, c1), (s2, c2)): ex = v
        if len(rows) % 20 == 0: print(f"  {label}: {len(rows)} pairs, {time.time() - t0:.0f} s", flush=True)
    d = pd.DataFrame(rows); d["matched_BH_FDR"] = bh(d.matched_p)
    print(f"  {label}: {len(d)} pairs in {time.time() - t0:.0f} s", flush=True)
    return d, ex

igc = {s: ig["data"][s]["comm"] for s in ig["data"]}
d_ig, ex = run(igc, "iGEA", example=(("Brass", "C4"), ("König", "C4")))
orc = {s: ora["res"][s]["comm"] for s in ora["res"]}
d_ora, _ = run(orc, "ORA")
# add the uniform null (s11f) for comparison
u = pd.read_csv(f"{H}/s11f_igea_pairs_adjusted.tsv", sep="\t")[["screen_1", "community_1", "screen_2", "community_2", "empirical_p", "BH_FDR"]]
d_ig = d_ig.merge(u.rename(columns={"empirical_p": "uniform_p", "BH_FDR": "uniform_BH_FDR"}), how="left",
                  on=["screen_1", "community_1", "screen_2", "community_2"]).sort_values("matched_p")
uo = pd.read_csv(f"{H}/s11f_ora_pairs_adjusted.tsv", sep="\t")[["screen_1", "community_1", "screen_2", "community_2", "empirical_p", "BH_FDR"]]
d_ora = d_ora.merge(uo.rename(columns={"empirical_p": "uniform_p", "BH_FDR": "uniform_BH_FDR"}), how="left",
                    on=["screen_1", "community_1", "screen_2", "community_2"]).sort_values("matched_p")
assert d_ig.uniform_p.notna().all() and len(d_ig) == 131 and len(d_ora) == 16
d_ig.to_csv(f"{H}/s11h_igea_pairs_matched{TAG}.tsv", sep="\t", index=False)
d_ora.to_csv(f"{H}/s11h_ora_pairs_matched{TAG}.tsv", sep="\t", index=False)
np.savez_compressed(f"{H}/s11h_example_null{TAG}.npz", brass_C4_konig_C4=ex)
summ = pd.DataFrame([dict(method=m, pairs=len(d), draws=N, matched_p_le_0_01=int((d.matched_p <= 0.01).sum()),
                          matched_BH_lt_0_05=int((d.matched_BH_FDR < 0.05).sum()), matched_BH_lt_0_10=int((d.matched_BH_FDR < 0.10).sum()),
                          uniform_BH_lt_0_05=int((d.uniform_BH_FDR < 0.05).sum())) for m, d in (("iGEA", d_ig), ("ORA", d_ora))])
summ.to_csv(f"{H}/s11h_summary{TAG}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(summ.to_string(index=False)); print(d_ig.head(15).to_string(index=False)); print(d_ora.head(4).to_string(index=False))
