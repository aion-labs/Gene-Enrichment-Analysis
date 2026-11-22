"""How many HIV ORA terms pass FDR under common variants of the BH denominator? (one-sided Fisher, 265 genes,
10 collections, sets restricted to the background, 10-600 genes). Reported terms always need overlap >= 3."""
import sys, os, logging
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s04_final_ora_vs_igea import P, R, BackgroundGeneSet, GeneConverter, GeneSet, GeneSetLibrary, FastCollection, exact_pvalue
bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
H = set(GeneSet(hiv, bg.genes, "HIV", hgcn=True, format=False).genes); B = sorted(bg.genes)
rows = []
for name, f in P:
    col = FastCollection.from_library(GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=name), B)
    x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in H if g in col.index]; n = x.sum(); k = col.A @ x
    p = np.array([exact_pvalue("Fisher's Exact Test", int(k[t]), int(col.K[t]), int(n), col.M) for t in range(len(col.term_names))])
    for t in range(len(p)): rows.append(dict(coll=name, k=int(k[t]), p=p[t]))
d = pd.DataFrame(rows)
def count(sub_mask, alpha, pooled):
    s = d[sub_mask].copy()
    if pooled: s["q"] = multipletests(s.p, method="fdr_bh")[1]
    else: s["q"] = s.groupby("coll").p.transform(lambda v: multipletests(v, method="fdr_bh")[1])
    return int(((s.q < alpha) & (s.k >= 3)).sum())
out = []
for label, mask in [("all sets 10-600 (paper; Enrichr-style)", d.k >= 0), ("sets with >= 1 input gene", d.k >= 1),
                    ("sets with >= 3 input genes (testable)", d.k >= 3)]:
    for pooled in (False, True):
        out.append({"BH denominator": label, "BH over": "all 10 collections pooled" if pooled else "each collection",
                    "m": int(mask.sum()), **{f"FDR < {a}": count(mask, a, pooled) for a in (0.05, 0.10, 0.25)}})
o = pd.DataFrame(out)
o.loc[len(o)] = {"BH denominator": "nominal p < 0.01 (no correction)", "BH over": "-", "m": int((d.k >= 3).sum()),
                 "FDR < 0.05": int(((d.p < 0.01) & (d.k >= 3)).sum()), "FDR < 0.10": None, "FDR < 0.25": None}
o.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "hiv_ora_fdr_variants.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250); print(o.to_string(index=False))
print(d[d.k >= 3].groupby("coll").size().to_dict())
