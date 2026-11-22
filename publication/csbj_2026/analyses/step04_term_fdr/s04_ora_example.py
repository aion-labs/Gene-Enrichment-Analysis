"""ORA vs empirical calibration, HIV x Reactome.
Real ORA: all Reactome sets (10-600) tested once on the 272 genes (Fisher, overlap >= 3), BH FDR over all tested sets.
Random lists: the same 10,000 lists as the per-run null (default_rng([2026, 272, i])), ORA on each.
 - empirical q (same logic as iGEA 'A'): mean # random ORA terms with p' <= p  /  # real ORA terms with p' <= p, monotone.
 - minP (same logic as iGEA 'B', single step): fraction of random lists whose BEST p <= p  (Westfall-Young, FWER).
Usage: python s04_ora_example.py <out_dir>"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging
from multiprocessing import get_context
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection, _fisher_two_sided_fast, exact_pvalue
W = {}
def init(col, B, n): W.update(col=col, B=B, n=n)
def ora_p(col, genes, exact=False):
    x = np.zeros(len(col.genes))
    for g in genes:
        i = col.index.get(g)
        if i is not None: x[i] = 1
    n = x.sum(); k = col.A @ x; cand = np.flatnonzero(k >= 3)
    if exact: p = np.array([exact_pvalue("Fisher's Exact Test", int(k[c]), int(col.K[c]), int(n), col.M) for c in cand])
    else: p = _fisher_two_sided_fast(k[cand], col.K[cand], n, col.M) if cand.size else np.array([])
    return cand, k, p, n
def one(i):
    g = [W["B"][j] for j in np.random.default_rng([2026, W["n"], i]).choice(len(W["B"]), W["n"], replace=False)]
    _, _, p, _ = ora_p(W["col"], g)
    return p[p < 0.01], (p.min() if p.size else 1.0)
if __name__ == "__main__":
    OUT = sys.argv[1]; N = 10000
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
    gs = GeneSet(hiv, bg.genes, "HIV", hgcn=False, format=False); B = sorted(set(bg.genes))
    col = FastCollection.from_library(GeneSetLibrary(R + "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt", name="Reactome"), B)
    # real ORA with exact p-values; BH over ALL tested sets (as in iGEA's Enrichment)
    x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in gs.genes if g in col.index]
    n = x.sum(); k = col.A @ x
    p_all = np.array([exact_pvalue("Fisher's Exact Test", int(k[c]), int(col.K[c]), int(n), col.M) for c in range(len(col.term_names))])
    fdr_all = multipletests(p_all, method="fdr_bh")[1]
    sel = np.flatnonzero((p_all < 0.01) & (k >= 3))
    real = pd.DataFrame({"Term": [col.term_names[c].replace("REACTOME_", "").replace("_", " ").lower() for c in sel],
                         "Overlap": k[sel].astype(int), "p": p_all[sel], "BH FDR": fdr_all[sel]}).sort_values("p").reset_index(drop=True)
    with get_context("spawn").Pool(os.cpu_count(), initializer=init, initargs=(col, B, gs.size)) as P:
        res = P.map(one, range(N), chunksize=100)
    rp = np.sort(np.concatenate([r[0] for r in res])); mins = np.array([r[1] for r in res])
    real["Real terms <= p"] = np.arange(1, len(real) + 1)
    real["Random terms <= p (per list)"] = np.searchsorted(rp, real.p.values, side="right") / N
    real["Empirical q"] = np.minimum.accumulate((real["Random terms <= p (per list)"] / real["Real terms <= p"]).clip(upper=1)[::-1])[::-1]
    real["minP (FWER)"] = [((mins <= p).sum() + 1) / (N + 1) for p in real.p]
    igea = {"regulation of glucokinase by glucokinase regulatory protein","asparagine n linked glycosylation","transcriptional regulation of white adipocyte differentiation","clathrin mediated endocytosis","nuclear rna decay","interleukin 12 signaling","retrograde transport at the trans golgi network","post translational modification synthesis of gpi anchored proteins","hiv transcription elongation"}
    real["iGEA term"] = real.Term.isin(igea).map({True: "yes", False: ""})
    real.to_csv(f"{OUT}/hiv_reactome_ora_calibration.tsv", sep="\t", index=False)
    print(f"ORA terms p<0.01 (overlap>=3): {len(real)}; random lists: {N}; random ORA terms p<0.01 per list: mean {len(rp)/N:.2f}, median {np.median([len(r[0]) for r in res]):.0f}")
    for lab, col_ in [("BH FDR", "BH FDR"), ("Empirical q", "Empirical q"), ("minP (FWER)", "minP (FWER)")]:
        print(f"  {lab:12} < 0.05: {(real[col_] < 0.05).sum():3d}   < 0.10: {(real[col_] < 0.10).sum():3d}")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 55)
    print(real.to_string(formatters={"p": "{:.1e}".format, "BH FDR": "{:.3f}".format, "Random terms <= p (per list)": "{:.3f}".format, "Empirical q": "{:.3f}".format, "minP (FWER)": "{:.4f}".format}))
