"""Data for Supplementary Table 2: gene sets per collection and HIV counts (current code, one universe, one-sided test).
Reads the HIV list from analyses/inputs and the ORA/iGEA counts from Supplementary Tables 4 and 5 (text).
Usage: python collection_table_data.py [<folder with Supp Table 04/05 .txt>] > collection_table_data.tsv"""
import pathlib as _pl, os
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
os.chdir(_REPO)
import sys, json, csv, statistics, logging; logging.disable(logging.WARNING)
sys.path.insert(0, "code")
from background_gene_set import BackgroundGeneSet
from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter
from fast_igea import FastCollection
SUPP = sys.argv[1] if len(sys.argv) > 1 else str(_PUB / "supplement")
alias = [a for a in json.load(open("data/libraries/alias.json")) if a["active"]]
B = sorted(BackgroundGeneSet("data/backgrounds/all_genes.txt").genes); Bs = set(B)
conv = GeneConverter()
hiv = [g for g in dict.fromkeys(m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(_PUB / "analyses" / "inputs" / "HIV_dependency_factors.symbols.txt") if l.strip()) if m) if g in Bs]
ora = list(csv.DictReader(open(SUPP + "/Supp Table 04 - ORA results.txt"), delimiter="\t"))
ig = list(csv.DictReader(open(SUPP + "/Supp Table 05 - iGEA results.txt"), delimiter="\t"))
rows = []
for a in alias:
    lib = GeneSetLibrary("data/libraries/" + a["file"], name=a["name"])
    col = FastCollection.from_library(lib, B)
    n_gmt = sum(1 for _ in open("data/libraries/" + a["file"]))
    sizes = [int(k) for k in col.K]
    x = set(hiv) & set(col.index)
    maxp = sum(1 for m in col.members if sum(1 for j in m if col.genes[j] in x) >= 3)
    o = [r for r in ora if r["Collection"] == a["name"]]
    i = [r for r in ig if r["Collection"] == a["name"]]
    rows.append(dict(Collection=a["name"], File=a["file"], **{
        "Gene sets in file": n_gmt, "Gene sets used (10-600 genes)": len(col.term_names),
        "Genes covered": col.M, "Median set size": int(statistics.median(sizes)),
        "HIV genes covered": len(x), "HIV: possible terms (>= 3 genes)": maxp,
        "HIV: ORA p < 0.01": sum(1 for r in o if float(r["p-value"]) < 0.01),
        "HIV: ORA FDR < 0.05": sum(1 for r in o if float(r["FDR"]) < 0.05),
        "HIV: iGEA selected": len(i), "HIV: iGEA kept (q < 0.10)": sum(1 for r in i if r["In network"] in ("True", "yes", "1"))}))
w = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
tot = {k: sum(r[k] for r in rows) for k in rows[0] if k not in ("Collection", "File", "Median set size", "Genes covered", "HIV genes covered")}
print("TOTAL", tot, "union genes:", len(set().union(*[set(FastCollection.from_library(GeneSetLibrary("data/libraries/" + a["file"], name=a["name"]), B).genes) for a in alias])), file=sys.stderr)
