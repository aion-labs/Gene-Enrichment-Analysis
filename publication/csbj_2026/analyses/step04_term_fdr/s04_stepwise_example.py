"""Step-wise (conditional) empirical p for each iteration, HIV x Reactome.
At iteration k the input has n_k remaining genes (the others were removed by earlier terms). Null for step k:
random lists of n_k genes drawn from the background minus the genes already removed, analysed exactly like the
real step (same collection, same M, min overlap 3); statistic = the BEST p-value among all gene sets.
p_step = (r + 1) / (N + 1), r = random lists whose best p <= observed p.
Usage: python s04_stepwise_example.py <out_dir> <N>"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging
from multiprocessing import get_context
import numpy as np, pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection
W = {}
def init(col, pool, n, seed): W.update(col=col, pool=pool, n=n, seed=seed)
def best(i):
    rng = np.random.default_rng([W["seed"], W["n"], i]); g = [W["pool"][j] for j in rng.choice(len(W["pool"]), W["n"], replace=False)]
    r = W["col"].iterate(g, p_threshold=1.0, min_overlap=3, max_iterations=1)
    return r[0][2] if r else 1.0
if __name__ == "__main__":
    OUT, N = sys.argv[1], int(sys.argv[2])
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
    gs = GeneSet(hiv, bg.genes, "HIV", hgcn=False, format=False); B = sorted(set(bg.genes))
    col = FastCollection.from_library(GeneSetLibrary(R + "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt", name="Reactome"), B)
    obs = col.iterate(gs.genes, 0.01, 3, None)
    q = pd.read_csv(f"{OUT}/hiv_reactome_qvalues.tsv", sep="\t").set_index("Iter")["q"]
    removed, remaining, rows = set(), set(gs.genes), []
    for it, t, p, rem in obs:
        pool = sorted(set(B) - removed); n = len(remaining)
        with get_context("spawn").Pool(os.cpu_count(), initializer=init, initargs=(col, pool, n, 2026 + it)) as P:
            bp = np.array(P.map(best, range(N), chunksize=100))
        r = int((bp <= p).sum())
        rows.append({"Iter": it, "Term": t.replace("REACTOME_", "").replace("_", " ").lower(), "Remaining genes n_k": n,
                     "Overlap": len(rem), "p": p, "Random best p <= p": r, "Step-wise empirical p": (r + 1) / (N + 1),
                     "Median random best p": float(np.median(bp)), "q (pooled, previous)": q.loc[it]})
        removed |= set(rem); remaining -= set(rem)
        print(rows[-1]["Iter"], rows[-1]["Term"][:45], n, f"{p:.2e}", r, f"{rows[-1]['Step-wise empirical p']:.4f}", flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/hiv_reactome_stepwise.tsv", sep="\t", index=False)
