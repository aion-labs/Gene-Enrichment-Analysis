"""Signal-adjusted empirical q (pi0-style), HIV, 10 collections.
Standard q assumes all 272 genes are noise (complete null). Here, genes assigned to terms already accepted are treated as
signal; chance terms are estimated from random lists of the remaining m0 = 272 - |signal genes| genes, drawn from the
background minus the signal genes. Iterate: accept terms with q_adj < alpha -> new signal set -> recompute, until stable.
Usage: python s04_q_signal_adjusted.py <out_dir> <alpha> <N>"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging
from multiprocessing import get_context
import numpy as np, pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection
P = [("Hallmark","h.all"),("Reactome","c2.cp.reactome"),("KEGG Legacy","c2.cp.kegg_legacy"),("KEGG MEDICUS","c2.cp.kegg_medicus"),
     ("BioCarta","c2.cp.biocarta"),("PID","c2.cp.pid"),("WikiPathways","c2.cp.wikipathways"),("GO BP","c5.go.bp"),("GO CC","c5.go.cc"),("GO MF","c5.go.mf")]
W = {}
def init(cols, pool, n): W.update(cols=cols, pool=pool, n=n)
def one(i):
    g = [W["pool"][j] for j in np.random.default_rng([2026, W["n"], i]).choice(len(W["pool"]), W["n"], replace=False)]
    return [(c, p) for c, col in enumerate(W["cols"]) for _, _, p, _ in col.iterate(g, 0.01, 3, None)]
def qvals(obs_p, rnd_p, N):
    op = np.array(obs_p); rp = np.sort(rnd_p)
    raw = np.array([np.searchsorted(rp, p, side="right") / N / (op <= p).sum() for p in op])
    o = np.argsort(op); q = np.empty(len(op)); q[o] = np.minimum.accumulate(np.clip(raw[o], 0, 1)[::-1])[::-1]; return q
if __name__ == "__main__":
    OUT, alpha, N = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
    gs = GeneSet(hiv, bg.genes, "HIV", hgcn=False, format=False); B = sorted(set(bg.genes))
    cols = [FastCollection.from_library(GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n), B) for n, f in P]
    real = [(c, t, p, set(g)) for c, col in enumerate(cols) for _, t, p, g in col.iterate(gs.genes, 0.01, 3, None)]
    signal = set(); history = []
    for rnd_i in range(4):
        m0 = gs.size - len(signal); pool = sorted(set(B) - signal)
        with get_context("spawn").Pool(os.cpu_count(), initializer=init, initargs=(cols, pool, m0)) as p:
            res = pd.DataFrame([x for r in p.map(one, range(N), chunksize=50) for x in r], columns=["c", "p"])
        q = np.empty(len(real))
        for c in range(len(cols)):
            idx = [i for i, r in enumerate(real) if r[0] == c]
            if idx: q[idx] = qvals([real[i][2] for i in idx], res.p[res.c == c].values, N)
        acc = [i for i in range(len(real)) if q[i] < alpha]
        new_signal = set().union(*[real[i][3] for i in acc]) if acc else set()
        history.append({"round": rnd_i, "m0 (noise genes)": m0, "chance terms/list": round(len(res) / N, 2), f"terms q<{alpha}": len(acc), "signal genes": len(new_signal)})
        print(history[-1], flush=True)
        if new_signal == signal: break
        signal = new_signal
    pd.DataFrame([{"Collection": P[r[0]][0], "Term": r[1], "p": r[2], "q_adjusted": q[i]} for i, r in enumerate(real)]).to_csv(f"{OUT}/hiv_q_signal_adjusted_alpha{alpha}.tsv", sep="\t", index=False)
    print("q_adj bands:", pd.cut(pd.Series(q), [-1, 0.05, 0.1, 0.25, 1.1]).value_counts().sort_index().to_dict())
