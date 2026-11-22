"""Two-layer statistics (q < 0.10, 10,000 random lists) for HIV with a subset of the 10 collections.
Random lists are analysed with the same subset, so the statistics belong to the network that is shown.
Usage: python subset_stats.py <subset: go|gobp|all10> <out_dir> [n_lists]"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection
import per_run_null as prn
ALL = [("H: Hallmark Gene Sets","h.all"),("C2: Reactome Pathways","c2.cp.reactome"),("C2: KEGG Legacy","c2.cp.kegg_legacy"),("C2: KEGG MEDICUS","c2.cp.kegg_medicus"),
       ("C2: BioCarta","c2.cp.biocarta"),("C2: Pathway Interaction Database","c2.cp.pid"),("C2: WikiPathways","c2.cp.wikipathways"),
       ("C5: Gene Ontology: Biological Process","c5.go.bp"),("C5: Gene Ontology: Cellular Component","c5.go.cc"),("C5: Gene Ontology: Molecular Function","c5.go.mf")]
SUB = {"all10": ALL, "go": ALL[7:], "gobp": ALL[7:8]}
if __name__ == "__main__":
    sub, out = sys.argv[1], sys.argv[2]; N = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
    gs = GeneSet(hiv, bg.genes, "HIV", hgcn=True, format=False)  # genes outside the background removed (one universe)
    libs = [GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n) for n, f in SUB[sub]]
    obs = [(c, t, it, p, g) for c, L in enumerate(libs) for it, t, p, g in FastCollection.from_library(L, bg.genes).iterate(gs.genes, 0.01, 3, 30)]
    table, null, terms, observed, attrs = prn.run_statistics(obs, libs, bg.genes, gs.size, n_lists=N, q_threshold=0.10, max_iterations=30,
                                                         cache_dir=os.path.expanduser('~/igea_cache'))
    params = dict(p_value_method="Fisher's Exact Test", p_threshold=0.01, min_overlap=3, min_term_size=10, max_term_size=600, max_iterations=30, background="all_genes")
    head = (f"Gene list: HIV dependency factors\nValidated input genes: {gs.size}\nCollections analysed (input and random lists): {len(libs)}\n"
            f"Largest cluster: {observed['cluster_size']} nodes, {observed['n_genes']} genes, {observed['n_terms']} terms, {observed['n_edges']} edges, "
            f"{observed['avg_edges_per_gene']:.4f} edges per gene, {observed['n_libraries']} collections\n\n")
    open(f"{out}/hiv_{sub}_statistical_report.txt", "w").write(head + prn.format_report(table, attrs, [L.name for L in libs], params, terms))
    terms.to_csv(f"{out}/hiv_{sub}_terms.tsv", sep="\t", index=False); null.to_csv(f"{out}/hiv_{sub}_random_lists.tsv", sep="\t", index=False)
    print(sub, attrs["seconds"], "s | kept", int(terms.Kept.sum()), "of", len(terms), "|", {k: observed[k] for k in ["n_genes", "n_terms", "n_edges"]})
    print(table[["Metric", "Observed", "Null Mean", "Z-score", "Random >= Observed", "Empirical p"]].round(4).to_string(index=False))
