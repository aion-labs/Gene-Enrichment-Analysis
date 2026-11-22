"""Supplementary Tables 4-6 (numbered 3-5 until 20 Nov) (HIV dependency factors, 10 collections), from the end-to-end CLI runs.
  Table 4: ORA (cli.py -m regular): every gene set tested once; p < 0.01, overlap >= 3; BH FDR per collection.
  Table 5: iGEA (cli.py -m iterative): terms by iteration with q-value, network membership and community, plus the
           same term's full-list overlap, p-value and BH FDR (as in ORA), for comparison.
  Table 6: statistical reports (10 collections, and the GO network of Figure 2).
Column header "Collection" (was "Library" before 17 Nov; values unchanged).
Usage: python build_supp_tables_4_6.py <ORA run folder> <step09_end_to_end> <figures/networks> <out_dir>"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, glob, logging
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection, exact_pvalue
import json
ORA, E9, NW, OUT = sys.argv[1:5]
ora = pd.read_csv(glob.glob(f"{ORA}/*/combined_regular_results.tsv")[0], sep="\t")
ora.rename(columns={"Library": "Collection"}).to_csv(f"{OUT}/Supp Table 04 - ORA results.txt", sep="\t", index=False)
run = glob.glob(f"{E9}/HIV_dependency_factors.symbols_*/")[0]
it = pd.read_csv(f"{run}/combined_iterative_results.tsv", sep="\t")
raw = pd.read_csv(glob.glob(f"{run}/*_term_qvalues.tsv")[0], sep="\t")      # raw term names, same order
# full-list statistics for every iGEA term (same computation as ORA)
alias = {a["name"]: a["file"] for a in json.load(open(R + "data/libraries/alias.json"))}
bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
H = set(GeneSet(hiv, bg.genes, "HIV", hgcn=True, format=False).genes); B = sorted(bg.genes)
full = {}
for lib in it.Library.unique():
    col = FastCollection.from_library(GeneSetLibrary(R + "data/libraries/" + alias[lib], name=lib), B)
    x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in H if g in col.index]; n = x.sum(); k = col.A @ x
    p = np.array([exact_pvalue("Fisher's Exact Test", int(k[t]), int(col.K[t]), int(n), col.M) for t in range(len(col.term_names))])
    fdr = multipletests(p, method="fdr_bh")[1]
    for t, name in enumerate(col.term_names):
        full[(lib, name)] = (f"{int(k[t])}/{int(col.K[t])}", p[t], fdr[t],
                             ", ".join(sorted(col.genes[j] for j in col.members[t] if col.genes[j] in H)))
key = list(zip(raw.Library, raw.Term))
assert len(key) == len(it) and all(a == b for a, b in zip(raw.Library, it.Library))
it["Full list overlap"] = [full[k][0] for k in key]
it["Full list p-value"] = [full[k][1] for k in key]
it["Full list FDR (BH)"] = [full[k][2] for k in key]
it["Full list overlapping genes"] = [full[k][3] for k in key]
it = it.rename(columns={"Overlap size": "Iteration overlap", "p-value": "Iteration p-value", "-log(p-value)": "Iteration -log10(p)",
                        "q-value": "q-value (empirical FDR)", "Genes": "Genes removed for next iteration"})
cols = ["Collection", "Iteration", "Term", "Description", "Iteration overlap", "Iteration p-value", "Iteration -log10(p)",
        "q-value (empirical FDR)", "In network", "Community", "Genes removed for next iteration",
        "Full list overlap", "Full list p-value", "Full list FDR (BH)", "Full list overlapping genes"]
it.rename(columns={"Library": "Collection"})[cols].to_csv(f"{OUT}/Supp Table 05 - iGEA results.txt", sep="\t", index=False)
hdr = lambda t: (t.replace("  Library                                   It ", "  Collection                                It ")  # newer CLI writes Collection
                 .replace("method=Fisher's Exact Test,", "method=Fisher's Exact Test (one-sided),"))
rep10 = hdr(open(glob.glob(f"{run}/*_statistical_report.txt")[0]).read())
repgo = hdr(open(f"{NW}/hiv_go_statistical_report.txt").read())
bar = "#" * 100
open(f"{OUT}/Supp Table 07 - Network statistical report.txt", "w").write(
    f"{bar}\nPART A. HIV dependency factors, 10 collections (Supplementary Figure 3; Results)\n{bar}\n\n{rep10}\n\n"
    f"{bar}\nPART B. HIV dependency factors, 3 Gene Ontology collections (Figure 2; random lists analysed with the same 3 collections)\n{bar}\n\n{repgo}\n")
print("Table 4:", len(ora), "terms;", int((ora.FDR < 0.05).sum()), "with FDR < 0.05")
print("Table 5:", len(it), "terms;", int(it["In network"].sum()), "in network; full-list FDR < 0.05:", int((it["Full list FDR (BH)"] < 0.05).sum()))
