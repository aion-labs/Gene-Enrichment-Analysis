"""Supp Fig 6 data: largest-cluster size and edges per gene of random lists of increasing size.
Same engine, collections, background and parameters as the CLI runs (one-sided Fisher, p < 0.01, overlap >= 3,
sets of 10-600 genes, <= 30 iterations). Usage: python s12a_size_trend.py <size> <n_lists> <q: none|0.10>
Appends rows to s12a_size_trend.tsv."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code")
logging.disable(logging.WARNING)
import pandas as pd
from pathlib import Path
from background_gene_set import BackgroundGeneSet
from gene_set_library import GeneSetLibrary
import cli
from per_run_null import generate_null
LIBS = ["h.all", "c2.cp.biocarta", "c2.cp.kegg_medicus", "c2.cp.pid", "c2.cp.reactome", "c2.cp.wikipathways",
        "c2.cp.kegg_legacy", "c5.go.bp", "c5.go.cc", "c5.go.mf"]
def main():
    size, n, q = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt")
    libs = []
    for f in LIBS:
        p = Path(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt")
        lib = GeneSetLibrary(str(p), name=cli.get_library_name_from_alias(p))
        lib.library = [t for t in lib.library if 10 <= t["size"] <= 600]
        lib.num_terms = len(lib.library); lib.unique_genes = lib.compute_unique_genes(); lib.size = len(lib.unique_genes)
        libs.append(lib)
    qt = None if q == "none" else float(q)
    null = generate_null(libs, bg.genes, size, n_lists=n, seed=2026, q_threshold=qt, max_iterations=30)
    null["list_size"] = size; null["q_filter"] = q
    out = Path(__file__).with_name("s12a_size_trend.tsv")
    null.to_csv(out, sep="\t", index=False, mode="a", header=not out.exists())
    print(size, q, len(null), null["cluster_size"].median(), null["avg_edges_per_gene"].median())
if __name__ == "__main__":
    main()
