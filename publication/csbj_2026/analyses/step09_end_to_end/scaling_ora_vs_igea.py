"""How do the ORA and iGEA networks grow as collections are added? (HIV, 265 genes)
ORA: terms at BH FDR < 0.05 (per collection), each term linked to all its genes.
iGEA: terms kept (q < 0.10) in the 10-collection run, each linked to its assigned genes.
Collections added in a fixed order: GO (3) -> +Reactome -> +KEGG Legacy, KEGG MEDICUS -> +WikiPathways, BioCarta,
PID -> +Hallmark. Prompt size = characters of the LLM prompt built with llm_export (same builder for both)."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import glob, os, sys
import pandas as pd
sys.path.insert(0, (str(_REPO) + "/code"))
from llm_export import semantic_dot, llm_prompt
from network_modules import analyze_network_modules
H = os.path.dirname(os.path.abspath(__file__))
ora_dir = sorted(glob.glob(f"{H}/ora_regular/HIV_dependency_factors.symbols_*"))[-1]
ig = pd.read_csv(sorted(glob.glob(f"{H}/HIV_dependency_factors.symbols_*/combined_iterative_results.tsv"))[-1], sep="\t")
GO = ["C5: Gene Ontology: Biological Process", "C5: Gene Ontology: Cellular Component", "C5: Gene Ontology: Molecular Function"]
STEPS = [("GO (3)", GO), ("+ Reactome (4)", ["C2: Reactome Pathways"]),
         ("+ KEGG (6)", ["C2: KEGG Legacy", "C2: KEGG MEDICUS"]),
         ("+ WikiPathways, BioCarta, PID (9)", ["C2: WikiPathways", "C2: BioCarta", "C2: Pathway Interaction Database"]),
         ("+ Hallmark (10)", ["H: Hallmark Gene Sets"])]
split = lambda s: [g.strip() for g in str(s).split(",") if g.strip()]
def ora_terms(lib):
    t = pd.read_csv(f"{ora_dir}/{lib}_regular_results.tsv", sep="\t"); t = t[t["FDR"] < 0.05]
    return [dict(library=lib, term=r.Term, p=r._6, q=r.FDR, genes=split(r.Genes)) for r in t.itertuples()]
def igea_terms(lib):
    t = ig[(ig["Library"] == lib) & (ig["In network"] == True)]
    return [dict(library=lib, term=r.Term, iteration=r.Iteration, p=r._6, q=r._8, genes=split(r.Genes)) for r in t.itertuples()]
def stats(terms):
    deg = {}
    for t in terms:
        for g in t["genes"]: deg[g] = deg.get(g, 0) + 1
    m = analyze_network_modules(terms)
    dot = semantic_dot(terms, m, None)
    return dict(terms=len(terms), genes=len(deg), edges=sum(deg.values()), max_degree=max(deg.values()),
                mean_degree=round(sum(deg.values()) / len(deg), 2), communities=m.summary.get("n_communities"),
                prompt_chars=len(llm_prompt(dot, terms, m, None)))
rows, libs = [], []
for name, add in STEPS:
    libs += add
    for meth, f in (("ORA", ora_terms), ("iGEA", igea_terms)):
        terms = [t for l in libs for t in f(l)]
        rows.append(dict(collections=name, method=meth, **stats(terms)))
o = pd.DataFrame(rows); o.to_csv(f"{H}/scaling_ora_vs_igea.tsv", sep="\t", index=False)
print(o.to_string(index=False))
