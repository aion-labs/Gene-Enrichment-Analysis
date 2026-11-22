"""Run the two-layer statistics (per_run_null.run_statistics) for a gene list with the 10 collections.
Usage: python run_two_layer.py <gene_list> <out_prefix> <n_lists> [q_threshold|none]"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, time, logging
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection
import per_run_null as prn
P = [("H: Hallmark Gene Sets","h.all"),("C2: Reactome Pathways","c2.cp.reactome"),("C2: KEGG Legacy","c2.cp.kegg_legacy"),("C2: KEGG MEDICUS","c2.cp.kegg_medicus"),
     ("C2: BioCarta","c2.cp.biocarta"),("C2: Pathway Interaction Database","c2.cp.pid"),("C2: WikiPathways","c2.cp.wikipathways"),
     ("C5: Gene Ontology: Biological Process","c5.go.bp"),("C5: Gene Ontology: Cellular Component","c5.go.cc"),("C5: Gene Ontology: Molecular Function","c5.go.mf")]
if __name__ == "__main__":
    gl, out, N = sys.argv[1], sys.argv[2], int(sys.argv[3]); qt = None if (len(sys.argv) > 4 and sys.argv[4] == "none") else float(sys.argv[4]) if len(sys.argv) > 4 else 0.10
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    genes = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(gl) if l.strip()) if m]
    gs = GeneSet(genes, bg.genes, "x", hgcn=False, format=False)
    libs = [GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n) for n, f in P]
    obs = [(c, t, it, p, g) for c, L in enumerate(libs) for it, t, p, g in FastCollection.from_library(L, bg.genes).iterate(gs.genes, 0.01, 3, None)]
    table, null, terms, observed, attrs = prn.run_statistics(obs, libs, bg.genes, gs.size, n_lists=N, q_threshold=qt)
    params = dict(p_value_method="Fisher's Exact Test", p_threshold=0.01, min_overlap=3, min_term_size=10, max_term_size=600, max_iterations=None, background="all_genes")
    open(out + "_report.txt", "w").write(prn.format_report(table, attrs, [n for n, _ in P], params, terms))
    terms.to_csv(out + "_terms.tsv", sep="\t", index=False); null.to_csv(out + "_random_lists.tsv", sep="\t", index=False)
    print({k: attrs[k] for k in ["seconds", "noise_genes", "expected_chance_terms", "grid_sizes"]})
    print("kept", int(terms.Kept.sum()), "of", len(terms), "| observed", {k: observed[k] for k in ["n_genes", "n_terms", "n_edges", "n_libraries"]})
    print("random lists: terms kept per list mean %.2f (all selected %.2f)" % (null.n_terms_kept.mean(), null.n_terms_total.mean()))
    print(table[["Metric", "Observed", "Null Mean", "Null SD", "Z-score", "Random >= Observed", "Empirical p"]].round(4).to_string(index=False))
