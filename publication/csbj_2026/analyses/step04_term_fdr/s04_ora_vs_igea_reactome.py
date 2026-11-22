"""HIV x Reactome: ORA terms with BH FDR < 0.05 vs iGEA terms with empirical q < 0.05.
For each ORA term: its overlapping HIV genes, how many are covered by the iGEA q<0.05 terms (genes assigned to them),
and the iGEA term (any q) that holds most of its genes."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging
import numpy as np, pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection
OUT = sys.argv[1]
bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
gs = GeneSet(hiv, bg.genes, "HIV", hgcn=False, format=False); B = sorted(set(bg.genes))
col = FastCollection.from_library(GeneSetLibrary(R + "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt", name="Reactome"), B)
nice = lambda t: t.replace("REACTOME_", "").replace("_", " ").lower()
ora = pd.read_csv(f"{OUT}/hiv_reactome_ora_calibration.tsv", sep="\t")
q = pd.read_csv(f"{OUT}/hiv_reactome_qvalues.tsv", sep="\t")
igea = col.iterate(gs.genes, 0.01, 3, None)
assigned = {nice(t): set(g) for _, t, _, g in igea}
qv = dict(zip(q.Term, q.q))
keep = {t for t in assigned if qv[t] < 0.05}
cov_keep = set().union(*[assigned[t] for t in keep])
hs = set(gs.genes); idx = {nice(n): i for i, n in enumerate(col.term_names)}
rows = []
for _, r in ora[ora["BH FDR"] < 0.05].iterrows():
    genes = {col.genes[j] for j in col.members[idx[r.Term]] if col.genes[j] in hs}
    best = max(assigned, key=lambda t: len(genes & assigned[t]))
    rows.append({"ORA term": r.Term, "ORA p": r.p, "BH FDR": r["BH FDR"], "HIV genes": len(genes),
                 "Covered by iGEA q<0.05 genes": len(genes & cov_keep),
                 "Main iGEA term (genes shared)": f"{best} ({len(genes & assigned[best])})" if genes & assigned[best] else "—",
                 "Status": "in both" if r.Term in keep else ("iGEA term, q>=0.05" if r.Term in assigned else "ORA only")})
d = pd.DataFrame(rows); d.to_csv(f"{OUT}/hiv_reactome_ora_fdr_vs_igea_q.tsv", sep="\t", index=False)
print("ORA FDR<0.05:", len(d), "| iGEA q<0.05:", len(keep), "| in both:", (d.Status == "in both").sum())
print("iGEA q<0.05 terms not ORA FDR<0.05:", sorted(keep - set(d["ORA term"])))
print("Genes in ORA FDR<0.05 terms:", len(set().union(*[{col.genes[j] for j in col.members[idx[t]] if col.genes[j] in hs} for t in d['ORA term']])),
      "| genes assigned to iGEA q<0.05 terms:", len(cov_keep), "| to all 9 iGEA terms:", len(set().union(*assigned.values())))
# group ORA-only terms by the iGEA term that holds most of their genes
d["main"] = d["Main iGEA term (genes shared)"].str.replace(r" \(\d+\)$", "", regex=True)
print(d.groupby("main").size().sort_values(ascending=False).to_string())
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 62)
print(d.drop(columns="main").to_string(index=False, formatters={"ORA p": "{:.1e}".format, "BH FDR": "{:.3f}".format}))
