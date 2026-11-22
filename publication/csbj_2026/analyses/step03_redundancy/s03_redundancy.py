"""Step 2 (R5-M3): redundancy between the 10 collections.
1) For every pair of collections: gene sets (size 10-600) in A that are identical to, or have Jaccard >= 0.9 with,
   a gene set in B.  2) Deduplicated GMTs: a set is dropped if it has J >= 0.9 with a set in a collection earlier
   in PRIORITY.  3) HIV network: cross-collection term pairs that are near-duplicates and the edges they carry.
Usage: python s03_redundancy.py <out_dir>"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, itertools, json
import numpy as np, pandas as pd
from scipy import sparse
R = (str(_REPO) + "/"); OUT = sys.argv[1]
PRIORITY = [("Hallmark","h.all"),("Reactome","c2.cp.reactome"),("KEGG Legacy","c2.cp.kegg_legacy"),("KEGG MEDICUS","c2.cp.kegg_medicus"),
            ("BioCarta","c2.cp.biocarta"),("PID","c2.cp.pid"),("WikiPathways","c2.cp.wikipathways"),("GO BP","c5.go.bp"),("GO CC","c5.go.cc"),("GO MF","c5.go.mf")]
import logging; logging.disable(logging.WARNING); sys.path.insert(0, R + "code")
from background_gene_set import BackgroundGeneSet
BG = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt").genes
def load(f):
    """Gene sets restricted to the background, then 10-600 genes (one universe, as in iGEA)."""
    out = []
    for line in open(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt"):
        p = line.rstrip("\n").split("\t")
        g = [x for x in dict.fromkeys(p[2:]) if x in BG]
        if 10 <= len(g) <= 600: out.append((p[0], p[1], g, frozenset(g)))
    return out
C = {n: load(f) for n, f in PRIORITY}
genes = sorted(set().union(*[s[3] for v in C.values() for s in v])); gi = {g: i for i, g in enumerate(genes)}
def mat(sets):
    r = np.repeat(np.arange(len(sets)), [len(s[3]) for s in sets]); c = [gi[g] for s in sets for g in s[3]]
    return sparse.csr_matrix((np.ones(len(c)), (r, c)), shape=(len(sets), len(genes)))
M = {n: mat(v) for n, v in C.items()}; size = {n: np.array([len(s[3]) for s in v]) for n, v in C.items()}
J90 = {}; ident = {}; best = {}
names = [n for n, _ in PRIORITY]
for a, b in itertools.permutations(names, 2):
    inter = (M[a] @ M[b].T).toarray()
    union = size[a][:, None] + size[b][None, :] - inter
    J = inter / union
    best[(a, b)] = J.max(1)
    J90[(a, b)] = J >= 0.9
    ident[(a, b)] = (J == 1.0)
rows = []
for a in names:
    r = {"Collection": a, "Gene sets (10-600)": len(C[a])}
    for b in names:
        r[b] = "" if a == b else int(J90[(a, b)].any(1).sum())
    r["Any other (J>=0.9)"] = int(np.any([J90[(a, b)].any(1) for b in names if b != a], axis=0).sum())
    r["Any other (identical)"] = int(np.any([ident[(a, b)].any(1) for b in names if b != a], axis=0).sum())
    rows.append(r)
tab = pd.DataFrame(rows); tab.to_csv(f"{OUT}/s03_redundancy_matrix_J90.csv", index=False)
print("Rows: collection A; columns: number of A's sets with Jaccard >= 0.9 to some set in collection B")
print(tab.to_string(index=False))
# deduplicated GMTs (drop a set if J >= 0.9 with a set in an EARLIER collection)
dropped = {}
for i, (a, f) in enumerate(PRIORITY):
    drop = np.zeros(len(C[a]), bool)
    for b, _ in PRIORITY[:i]: drop |= J90[(a, b)].any(1)
    dropped[a] = int(drop.sum())
    with open(f"{OUT}/dedup_gmt/{f}.dedup.gmt", "w") as fh:
        for s, d in zip(C[a], drop):
            if not d: fh.write("\t".join([s[0], s[1]] + s[2]) + "\n")
print("\nDeduplicated collections — sets removed:", dropped, "total", sum(dropped.values()))
json.dump({"removed_by_dedup": dropped, "priority": names}, open(f"{OUT}/s03_dedup_summary.json", "w"), indent=2)
# HIV network: near-duplicate term pairs across collections
import glob
it = pd.read_csv(glob.glob((str(_PUB) + "/analyses/step09_end_to_end/HIV_dependency_factors.symbols_*/HIV_dependency_factors.symbols_term_qvalues.tsv"))[0], sep="\t")
it = it[it.Kept.astype(bool)]  # final HIV network (q < 0.10)
lib2 = {"H: Hallmark Gene Sets":"Hallmark","C2: Reactome Pathways":"Reactome","C2: KEGG Legacy":"KEGG Legacy","C2: KEGG MEDICUS":"KEGG MEDICUS","C2: BioCarta":"BioCarta",
        "C2: Pathway Interaction Database":"PID","C2: WikiPathways":"WikiPathways","C5: Gene Ontology: Biological Process":"GO BP","C5: Gene Ontology: Cellular Component":"GO CC","C5: Gene Ontology: Molecular Function":"GO MF"}
setof = {s[0]: (n, s[3]) for n, v in C.items() for s in v}
raw = {}
for _, r in it.iterrows():
    t = r.Term
    cand = [k for k in setof if k == t]
    raw[r.Term] = (lib2[r.Library], cand[0] if cand else None, set(str(r.Genes).split(", ")))
pairs = []
for (t1, (c1, k1, g1)), (t2, (c2, k2, g2)) in itertools.combinations(raw.items(), 2):
    if c1 == c2 or not k1 or not k2: continue
    s1, s2 = setof[k1][1], setof[k2][1]; j = len(s1 & s2) / len(s1 | s2)
    if j >= 0.5: pairs.append({"Term 1": t1, "Coll 1": c1, "Term 2": t2, "Coll 2": c2, "Jaccard (gene sets)": round(j, 3), "Shared network genes": len(g1 & g2)})
pp = pd.DataFrame(pairs).sort_values("Jaccard (gene sets)", ascending=False) if pairs else pd.DataFrame()
pp.to_csv(f"{OUT}/s03_hiv_similar_term_pairs.csv", index=False)
unmatched = [t for t, v in raw.items() if not v[1]]
print(f"\nHIV network: {len(raw)} terms ({len(unmatched)} not matched to GMT names: {unmatched[:3]})")
print("Cross-collection term pairs with gene-set Jaccard >= 0.5:"); print(pp.to_string(index=False) if len(pp) else "none")
