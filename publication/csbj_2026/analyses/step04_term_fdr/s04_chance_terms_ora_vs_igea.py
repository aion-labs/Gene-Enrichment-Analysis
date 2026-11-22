"""Chance terms per random list (N_GENES = validated HIV input in the background, 265; same lists as the per-run null): ORA (p<0.01, overlap>=3) vs iGEA, per collection."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging
from multiprocessing import get_context
import numpy as np, pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set_library import GeneSetLibrary
from fast_igea import FastCollection, _fisher_one_sided_fast
P = [("Hallmark","h.all"),("Reactome","c2.cp.reactome"),("KEGG Legacy","c2.cp.kegg_legacy"),("KEGG MEDICUS","c2.cp.kegg_medicus"),
     ("BioCarta","c2.cp.biocarta"),("PID","c2.cp.pid"),("WikiPathways","c2.cp.wikipathways"),("GO BP","c5.go.bp"),("GO CC","c5.go.cc"),("GO MF","c5.go.mf")]
W = {}
def init(cols, B, n): W.update(cols=cols, B=B, n=n)
def one(i):
    g = [W["B"][j] for j in np.random.default_rng([2026, W["n"], i]).choice(len(W["B"]), W["n"], replace=False)]
    out = []
    for col in W["cols"]:
        x = np.zeros(len(col.genes))
        for gg in g:
            j = col.index.get(gg)
            if j is not None: x[j] = 1
        k = col.A @ x; c = np.flatnonzero(k >= 3)
        n_ora = int((_fisher_one_sided_fast(k[c], col.K[c], x.sum(), col.M) < 0.01).sum()) if c.size else 0
        out.append((n_ora, len(col.iterate(g, 0.01, 3, 30))))
    return out
if __name__ == "__main__":
    OUT, N = sys.argv[1], int(sys.argv[2])
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); B = sorted(set(bg.genes))
    cols = [FastCollection.from_library(GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n), B) for n, f in P]
    with get_context("spawn").Pool(os.cpu_count(), initializer=init, initargs=(cols, B, 265)) as p:
        res = np.array(p.map(one, range(N), chunksize=50))           # N x 10 x 2
    d = pd.DataFrame({"Collection": [n for n, _ in P], "ORA chance terms / list": res[:, :, 0].mean(0).round(2),
                      "iGEA chance terms / list": res[:, :, 1].mean(0).round(2)})
    d.loc[len(d)] = ["Total", res[:, :, 0].sum(1).mean().round(1), res[:, :, 1].sum(1).mean().round(1)]
    d.to_csv(f"{OUT}/chance_terms_ora_vs_igea.tsv", sep="\t", index=False); print(d.to_string(index=False))
    print("ORA total per list: median", np.median(res[:, :, 0].sum(1)), "95th", np.quantile(res[:, :, 0].sum(1), .95),
          "| iGEA total: median", np.median(res[:, :, 1].sum(1)), "95th", np.quantile(res[:, :, 1].sum(1), .95))
