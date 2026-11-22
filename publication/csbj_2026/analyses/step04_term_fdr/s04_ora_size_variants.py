"""HIV ORA (one-sided Fisher, BH within each collection over all sets in the size range) with and without the
10-600 gene-set size filter; terms reported at FDR < 0.05 with overlap >= 3 or >= 1."""
import sys, os
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s04_final_ora_vs_igea import P, R, BackgroundGeneSet, GeneConverter, GeneSet, GeneSetLibrary, FastCollection, exact_pvalue
bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
H = set(GeneSet(hiv, bg.genes, "HIV", hgcn=True, format=False).genes); B = sorted(bg.genes)
rows = []
for lo, hi, lab in [(10, 600, "10-600 (paper)"), (1, 10**6, "no size filter"), (5, 10**6, ">= 5 genes"), (15, 500, "15-500 (GSEA default)")]:
    tot = dict(sizes=lab, sets=0, fdr05_k3=0, fdr05_k1=0, genes_k3=set())
    for name, f in P:
        col = FastCollection(GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=name).library, B, name, lo, hi)
        x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in H if g in col.index]; n = x.sum(); k = col.A @ x
        p = np.array([exact_pvalue("Fisher's Exact Test", int(k[t]), int(col.K[t]), int(n), col.M) for t in range(len(col.term_names))])
        q = multipletests(p, method="fdr_bh")[1]
        tot["sets"] += len(p); tot["fdr05_k3"] += int(((q < 0.05) & (k >= 3)).sum()); tot["fdr05_k1"] += int((q < 0.05).sum())
        for t in np.where((q < 0.05) & (k >= 3))[0]:
            tot["genes_k3"] |= {col.genes[j] for j in col.members[t] if col.genes[j] in H}
    tot["genes_k3"] = len(tot["genes_k3"]); rows.append(tot)
o = pd.DataFrame(rows); o.columns = ["Gene-set sizes", "Sets tested", "FDR<0.05, overlap>=3", "FDR<0.05, any overlap", "Genes (overlap>=3)"]
o.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "hiv_ora_size_variants.tsv"), sep="\t", index=False)
print(o.to_string(index=False))
