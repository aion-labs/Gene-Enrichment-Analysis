"""Step 2 sensitivity (two-level statistics, q < 0.10, as in the final method): HIV network and per-run null (a) with deduplicated collections, (b) leaving out
one collection at a time. Observed networks use the fast engine (identical to IterativeEnrichment).
Usage: python s03b_sensitivity.py <out_dir> <mode: dedup|loo> <n_lists>"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, json, logging
import pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter
from fast_igea import FastCollection
import per_run_null as prn
OUT, MODE, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
P = [("Hallmark","h.all"),("Reactome","c2.cp.reactome"),("KEGG Legacy","c2.cp.kegg_legacy"),("KEGG MEDICUS","c2.cp.kegg_medicus"),
     ("BioCarta","c2.cp.biocarta"),("PID","c2.cp.pid"),("WikiPathways","c2.cp.wikipathways"),("GO BP","c5.go.bp"),("GO CC","c5.go.cc"),("GO MF","c5.go.mf")]
if __name__ == "__main__":
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
    gs = GeneSet(hiv, bg.genes, "HIV", hgcn=True, format=False)  # one universe
    def lib(n, f, dedup=False):
        path = f"{OUT}/dedup_gmt/{f}.dedup.gmt" if dedup else R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt"
        return GeneSetLibrary(path, name=n)
    configs = {"all 10 (deduplicated)": [lib(n, f, True) for n, f in P]} if MODE == "dedup" else \
              {f"without {x}": [lib(n, f) for n, f in P if n != x] for x, _ in P}
    rows = []
    for name, libs in configs.items():
        rec = [(c, t, it, p, g) for c, L in enumerate(libs) for it, t, p, g in FastCollection.from_library(L, bg.genes).iterate(gs.genes, 0.01, 3, 30)]
        table, null, terms, obs, attrs = prn.run_statistics(rec, libs, bg.genes, gs.size, n_lists=N, seed=2026, q_threshold=0.10,
                                                            max_iterations=30, cache_dir=os.path.expanduser("~/igea_cache"))
        t = table.set_index("Metric")
        rows.append({"Configuration": name, "iGEA terms": len(rec), "Kept (q < 0.10)": obs["n_terms_kept"], "Cluster genes": obs["n_genes"], "Cluster terms": obs["n_terms"],
                     "Cluster edges": obs["n_edges"], "Edges/gene": round(obs["avg_edges_per_gene"], 3), "Collections in cluster": obs["n_libraries"],
                     **{f"z {m}": round(t.loc[m, "Z-score"], 2) for m in t.index},
                     **{f"p {m}": round(t.loc[m, "Empirical p"], 5) for m in t.index}, "Random lists": N})
        print(name, rows[-1]["Cluster genes"], rows[-1]["Cluster terms"], rows[-1]["Cluster edges"], {k: v for k, v in rows[-1].items() if k.startswith("p ")}, flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/s03_sensitivity_{MODE}.csv", index=False)
