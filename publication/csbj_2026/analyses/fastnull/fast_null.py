"""Generate a random-list null with the fast engine and summarise the largest cluster per list.
Usage: python fast_null.py <size> <n_lists> <seed> <out_csv> [background_file]
Sampling: uniform without replacement from the background (default all_genes.txt, as the stored null).
Seeded with numpy Generator on a SORTED gene list, so lists are exactly reproducible."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, time
from multiprocessing import Pool
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fast_igea import Collection
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code")
import logging; logging.disable(logging.WARNING)
from network_connectivity_benchmark import NetworkConnectivityAnalyzer
F = {"H: Hallmark Gene Sets":"h.all","C2: CP: BioCarta":"c2.cp.biocarta","C2: CP: KEGG MEDICUS":"c2.cp.kegg_medicus",
     "C2: CP: Pathway Interaction Database":"c2.cp.pid","C2: CP: Reactome Pathways":"c2.cp.reactome","C2: CP: WikiPathways":"c2.cp.wikipathways",
     "C2: CP: KEGG Legacy":"c2.cp.kegg_legacy","C5: GO: Biological Process":"c5.go.bp","C5: GO: Cellular Component":"c5.go.cc","C5: GO: Molecular Function":"c5.go.mf"}
size, nl, seed, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
bgf = sys.argv[5] if len(sys.argv) > 5 else R + "data/backgrounds/all_genes.txt"
BG = sorted({l.strip().upper() for l in open(bgf) if l.strip()})
COLS = None
def init():
    global COLS; COLS = {lib: Collection(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", BG) for lib, f in F.items()}
def one(i):
    rng = np.random.default_rng([seed, size, i]); g = [BG[j] for j in rng.choice(len(BG), size, replace=False)]
    res = []
    for lib, c in COLS.items():
        for it, term, p, removed in c.igea(g):
            res.append({"Term": term, "Iteration": it, "Library": lib, "Genes removed for next iteration": removed})
    a = NetworkConnectivityAnalyzer(); a.add_igea_results(res); a.compute_metrics(include_library_diversity=True); cl = a.get_clusters()
    c0 = cl[0] if cl else {"size": 0, "n_genes": 0, "n_terms": 0, "n_edges": 0, "avg_edges_per_gene": 0}
    return {"list": i, "size": c0["size"], "genes": c0["n_genes"], "terms": c0["n_terms"], "edges": c0["n_edges"],
            "edges_per_gene": c0.get("avg_edges_per_gene", 0), "n_terms_total": len(res)}
if __name__ == "__main__":
    t = time.time()
    NJ = int(os.environ.get("N_JOBS", os.cpu_count()))
    with Pool(NJ, initializer=init) as p: rows = p.map(one, range(nl), chunksize=10)
    pd.DataFrame(rows).to_csv(out, index=False); print(f"{nl} lists of {size} in {time.time()-t:.1f}s on {NJ} cores")
