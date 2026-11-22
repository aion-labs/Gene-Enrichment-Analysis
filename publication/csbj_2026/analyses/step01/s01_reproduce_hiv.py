"""Step 1 — reproduce the published HIV numbers with the frozen iGEA code (read-only use of the repo).
Run: python s01_reproduce_hiv.py <repo_path> <out_dir>
Parameters = the values used for the paper / permutation data: Fisher's exact test, p < 0.01,
min overlap 3, term size 10-600, background all_genes.txt, max iterations 30 (benchmark cap)."""
import sys, json, logging, time
from pathlib import Path
REPO = Path(sys.argv[1]).resolve(); OUT = Path(sys.argv[2]).resolve(); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "code"))
logging.basicConfig(level=logging.WARNING)
for n in ["enrichment","iter_enrichment","gene_set","gene_set_library","background_gene_set","gene_converter","parallel_null_distribution","network_connectivity_benchmark"]:
    logging.getLogger(n).setLevel(logging.WARNING)
import pandas as pd
from gene_converter import GeneConverter
from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from enrichment import Enrichment
from iter_enrichment import IterativeEnrichment
from network_connectivity_benchmark import NetworkConnectivityAnalyzer
from parallel_null_distribution import compute_null_distribution_from_raw_permutations

LIBS = {  # display name in paper tables : (gmt file, name in permutation data)
 "H: Hallmark Gene Sets": ("h.all.v2025.1.Hs.symbols.gmt", "H: Hallmark Gene Sets"),
 "C2: BioCarta": ("c2.cp.biocarta.v2025.1.Hs.symbols.gmt", "C2: CP: BioCarta"),
 "C2: KEGG MEDICUS": ("c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt", "C2: CP: KEGG MEDICUS"),
 "C2: Pathway Interaction Database": ("c2.cp.pid.v2025.1.Hs.symbols.gmt", "C2: CP: Pathway Interaction Database"),
 "C2: Reactome Pathways": ("c2.cp.reactome.v2025.1.Hs.symbols.gmt", "C2: CP: Reactome Pathways"),
 "C2: WikiPathways": ("c2.cp.wikipathways.v2025.1.Hs.symbols.gmt", "C2: CP: WikiPathways"),
 "C2: KEGG Legacy": ("c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt", "C2: CP: KEGG Legacy"),
 "C5: Gene Ontology: Biological Process": ("c5.go.bp.v2025.1.Hs.symbols.gmt", "C5: GO: Biological Process"),
 "C5: Gene Ontology: Cellular Component": ("c5.go.cc.v2025.1.Hs.symbols.gmt", "C5: GO: Cellular Component"),
 "C5: Gene Ontology: Molecular Function": ("c5.go.mf.v2025.1.Hs.symbols.gmt", "C5: GO: Molecular Function"),
}
P, MINOV, MINS, MAXS, MAXIT, METHOD = 0.01, 3, 10, 600, 30, "Fisher's Exact Test"
t0 = time.time()
assert Path.cwd().joinpath("data/recent_release/Homo_sapiens.gene_info").exists(), "run from a working dir containing data/recent_release (decompressed NCBI files)"
bg = BackgroundGeneSet(str(REPO / "data/backgrounds/all_genes.txt"))
conv = GeneConverter()
raw = [l.strip() for l in open(REPO / "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()]
mapped = [m for m in (conv.validate_and_map_symbol(g) for g in raw) if m]
gs = GeneSet(mapped, bg.genes, "HIV_dependency_factors.symbols", hgcn=False, format=False)
assert gs.size > 0, "no input genes mapped"
info = {"input_lines": len(raw), "mapped": len(mapped), "gene_set_size": gs.size, "background_size": bg.size}

ora_rows, it_rows = [], []
for disp, (gmt, _) in LIBS.items():
    lib = GeneSetLibrary(str(REPO / "data/libraries" / gmt), name=disp)
    enr = Enrichment(gs, lib, bg, MINS, MAXS, METHOD, use_multiprocessing=False)
    for r in enr.results:
        ov = int(r["overlap_size"].split("/")[0])
        ora_rows.append({"Library": disp, "Term": r["term"], "overlap": ov, "p": r["p-value"], "fdr": r["fdr"]})
    ie = IterativeEnrichment(gs, lib, bg, MINS, MAXS, METHOD, p_threshold=P, max_iterations=None,
                             min_overlap=MINOV, use_multiprocessing=False, save_iteration_results=False)
    for r in ie.results:
        it_rows.append({k: r[k] for k in ["Library","Iteration","Term","iteration p-value","iteration overlapping genes"]}
                       | {"Genes removed for next iteration": ", ".join(r["Genes removed for next iteration"])})
ora = pd.DataFrame(ora_rows); it = pd.DataFrame(it_rows)
ora.to_csv(OUT / "s01_ora_all_tested.tsv", sep="\t", index=False); it.to_csv(OUT / "s01_igea_terms.tsv", sep="\t", index=False)
ora_sig = ora[(ora.p < P) & (ora.overlap >= MINOV)]
info |= {"ora_terms_p<0.01_overlap>=3": len(ora_sig), "ora_terms_fdr<0.05": int((ora_sig.fdr < 0.05).sum()),
         "ora_fdr<0.05_any_overlap": int((ora.fdr < 0.05).sum()),
         "igea_terms": len(it), "igea_max_iteration": int(it.Iteration.max()),
         "igea_terms_iter<=30": int((it.Iteration <= MAXIT).sum())}
# network (same logic as cli.compute_benchmark_for_cli)
an = NetworkConnectivityAnalyzer()
sel = it[(it.Iteration <= MAXIT) & (it["iteration p-value"] <= P)]
an.add_igea_results([{"Term": r.Term, "Iteration": r.Iteration, "Library": r.Library,
    "Genes removed for next iteration": [g for g in r["Genes removed for next iteration"].split(", ") if g]} for _, r in sel.iterrows()])
an.compute_metrics(include_library_diversity=True)
lc = an.get_clusters()[0]
info |= {"largest_cluster": {k: lc[k] for k in ["size","n_genes","n_terms","n_edges"]} | {"avg_edges_per_gene": round(lc.get("avg_edges_per_gene", 0), 4), "n_libraries": lc.get("n_libraries")}}
# null (same function the CLI/app uses)
null = compute_null_distribution_from_raw_permutations(REPO / "permutations/merged_permutation_results.tsv", gs.size,
        [v[1] for v in LIBS.values()], P, MAXIT)
ns, size_key = null
obs = {"largest_component_size": lc["size"], "largest_cluster_genes": lc["n_genes"], "largest_cluster_terms": lc["n_terms"],
       "largest_cluster_edges": lc["n_edges"], "largest_cluster_avg_edges_per_gene": lc.get("avg_edges_per_gene", 0)}
info["null_size_key"] = size_key
info["z"] = {m: {"obs": v, "null_mean": round(ns[m]["mean"], 4), "null_sd": round(ns[m]["std"], 4), "n": ns[m]["n"],
                 "z": round((v - ns[m]["mean"]) / ns[m]["std"], 2)} for m, v in obs.items() if m in ns}
info["runtime_s"] = round(time.time() - t0, 1)
json.dump(info, open(OUT / "s01_reproduction.json", "w"), indent=2, default=str)
print(json.dumps(info, indent=2, default=str))
