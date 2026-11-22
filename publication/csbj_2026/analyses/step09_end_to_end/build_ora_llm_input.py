"""ORA baseline for the LLM comparison (R2-4): the same prompt as iGEA, applied to the FDR-filtered ORA network.

Input: CLI regular (ORA) run of the HIV list (ora_regular/, one-sided Fisher, 265 genes), three GO collections,
terms with Benjamini-Hochberg FDR < 0.05 (within collection) = the ORA side of Figure 1 / Figure 2.
Network: each term linked to ALL of its overlapping input genes (no iterative removal). Communities: the same
Leiden procedure as for iGEA. Prompt: llm_export.llm_prompt with only the method-specific sentences changed
(introduction, statistics, gene degree, q attribute); every replacement is asserted.
Usage: build_ora_llm_input.py [go|all]. Output folder: ora_go_fdr0.05/ or ora_all_fdr0.05/ (network_llm.dot, llm_prompt.txt, network_modules_*.tsv, ora_terms.tsv).
"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import glob, os, sys
import pandas as pd
sys.path.insert(0, (str(_REPO) + "/code"))
from llm_export import semantic_dot, llm_prompt, INTRO
from network_modules import analyze_network_modules, write_network_modules

HERE = os.path.dirname(os.path.abspath(__file__))
run = sorted(glob.glob(os.path.join(HERE, "ora_regular", "HIV_dependency_factors.symbols_*")))[-1]
GO = ["C5: Gene Ontology: Biological Process", "C5: Gene Ontology: Cellular Component",
      "C5: Gene Ontology: Molecular Function"]
ALL = ["H: Hallmark Gene Sets", "C2: BioCarta", "C2: KEGG MEDICUS", "C2: Pathway Interaction Database",
       "C2: Reactome Pathways", "C2: WikiPathways", "C2: KEGG Legacy"] + GO
SET = sys.argv[1] if len(sys.argv) > 1 else "go"          # go (Figure 2 comparison) | all (10 collections)
LIBS = GO if SET == "go" else ALL
rows = []
for lib in LIBS:
    t = pd.read_csv(os.path.join(run, f"{lib}_regular_results.tsv"), sep="\t")
    t = t[t["FDR"] < 0.05]
    for r in t.itertuples():
        rows.append(dict(library=lib, term=r.Term, p=r._6, q=r.FDR,
                         genes=[g.strip() for g in str(r.Genes).split(",") if g.strip()]))
out = os.path.join(HERE, "ora_go_fdr0.05" if SET == "go" else "ora_all_fdr0.05"); os.makedirs(out, exist_ok=True)
pd.DataFrame([dict(Collection=r["library"], Term=r["term"], p=r["p"], FDR=r["q"], Genes=", ".join(r["genes"]))
              for r in rows]).to_csv(os.path.join(out, "ora_terms.tsv"), sep="\t", index=False)
mods = analyze_network_modules(rows)
write_network_modules(mods, out)
dot = semantic_dot(rows, mods, None)
prompt = llm_prompt(dot, rows, mods, None)
n_genes = len({g for r in rows for g in r["genes"]})

ORA_INTRO = ("Below is a gene set enrichment network computed from a gene list by over-representation analysis (ORA): "
             "every gene set (= term) of each collection was tested against the full gene list, and each significant "
             "term is linked to all of its genes in the list. Thus a gene can be linked to several terms of the same "
             "collection. The network of genes and gene sets is given in DOT format, preceded by its statistics and "
             "communities.")
STATS_OLD = prompt[prompt.index("**STATISTICAL CONTEXT:**"):prompt.index("**COMMUNITIES:**")]
STATS_NEW = ("**STATISTICAL CONTEXT:**\nThe input contained 265 validated genes. Each term was tested with a one-sided "
             "Fisher's exact test against the full gene list; the "
             f"{len(rows)} terms with a Benjamini-Hochberg false discovery rate below 0.05 (within their collection) "
             f"form this network ({n_genes} genes). No network-level test was computed for this network. Do not use "
             "p-values or FDR values to rank biological importance.\n\n")
REPL = [
    (INTRO, ORA_INTRO),
    (STATS_OLD, STATS_NEW),
    ("Gene degree is the number of terms a gene connects to; in iGEA this equals the number of collections that "
     "place the gene in the network.", "Gene degree is the number of terms a gene connects to."),
    ("(number of terms the gene connects to; in iGEA, the number of collections that place it in the network)",
     "(number of terms the gene connects to)"),
    ("`iteration`, `p` (enrichment p-value in its iteration), `q` (term q-value)",
     "`p` (enrichment p-value) and `q` (Benjamini-Hochberg false discovery rate)"),
    ("Start with whether the network is distinguishable from random gene lists, then the overall biological picture",
     "Start with the overall biological picture"),
]
for old, new in REPL:
    assert prompt.count(old) == 1, old[:80]
    prompt = prompt.replace(old, new)
assert "iGEA" not in prompt and "iteration" not in prompt.lower().replace("iterative", ""), "iGEA wording left"
open(os.path.join(out, "network_llm.dot"), "w").write(dot)
open(os.path.join(out, "llm_prompt.txt"), "w").write(prompt)
print(f"{len(rows)} terms, {n_genes} genes, {mods.summary.get('n_communities')} communities "
      f"(modularity {mods.summary.get('modularity')}); prompt {len(prompt):,} chars -> {out}")
