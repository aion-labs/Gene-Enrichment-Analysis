"""Supplementary Table 6 (tab-separated text): ORA and iGEA networks as gene-set collections are added (HIV, hypoxia, vaccination).
Input: analyses/step09_end_to_end/scaling_all.tsv (written by scaling_all.py). HTML: build_supp_tables_04_07_html.py.
Usage: python build_supp_table_06_txt.py [<scaling_all.tsv>] [<out dir>]"""
import sys, os
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "analyses", "step09_end_to_end", "scaling_all.tsv")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "..", "supplement")
HEADER = 'Supplementary Table 6. ORA and iGEA networks as gene-set collections are added\n\nCollections are added in a fixed order: the three Gene Ontology collections, then Reactome, KEGG Legacy and KEGG MEDICUS,\nWikiPathways, BioCarta and Pathway Interaction Database, and Hallmark (10 collections).\nORA: hypergeometric test (upper tail) on the full gene list; gene sets restricted to the background and to 10-600 genes (the same\nsets as iGEA); Benjamini-Hochberg FDR within each collection; terms with FDR < 0.05 and at least 3 input genes, each linked to all\nits input genes. "ORA (all set sizes)": the same without the 10-600 filter, as ORA is usually run.\niGEA: terms kept (empirical q < 0.10, 10,000 random gene lists) in the 10-collection analysis, restricted to the collections of\neach row; each gene is assigned to at most one term per collection.\n(For HIV with GO alone, the separately tested GO network of Figure 2 and Figure 3B has 26 terms and\n138 genes; the GO rows below are the GO terms of the 10-collection analysis, 27 terms and 144 genes.)\nEdges: gene-term links. Max gene degree: largest number of terms linked to one gene. Communities and modularity: Leiden\n(igraph, modularity) on the term network, terms linked by shared genes, best of 100 seeds (the procedure of the iGEA software);\nmodularity describes how cleanly the terms separate into communities (descriptive, not a test).\n'
STEPS = ["GO (3)", "+ Reactome (4)", "+ KEGG Legacy, KEGG MEDICUS (6)", "+ WikiPathways, BioCarta, PID (9)", "+ Hallmark (10)", "All 10, all set sizes"]
d = pd.read_csv(SRC, sep="\t").drop_duplicates()
d["_ds"] = d.dataset.map({"HIV": 0, "Hypoxia": 1, "Vaccination": 2}); d["_st"] = d.collections.map(STEPS.index)
d["_m"] = d.method.map({"ORA": 0, "iGEA": 1, "ORA (all set sizes)": 2})
d = d.sort_values(["_ds", "_st", "_m"])
rows = ["\t".join(["Dataset", "Input genes", "Collections", "Method", "Terms", "Genes", "Edges", "Edges per gene", "Max gene degree",
                   "Communities", "Modularity"])]
for r in d.itertuples():
    rows.append("\t".join(str(v) for v in [r.dataset, r.input_genes, r.collections, r.method, r.terms, r.genes, r.edges,
                                           round(r.edges / r.genes, 2), r.max_gene_degree, r.communities, round(r.modularity, 2)]))
open(os.path.join(OUT, "Supp Table 06 - ORA vs iGEA as collections are added.txt"), "w").write(HEADER + "\n" + "\n".join(rows) + "\n")
print("written", len(rows) - 1, "rows")
