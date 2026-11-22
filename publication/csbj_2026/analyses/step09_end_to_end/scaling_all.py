"""ORA vs iGEA networks as collections are added, for HIV, hypoxia and vaccination (Supplementary Table 14).
ORA: one-sided Fisher on the full list, gene sets restricted to the background, 10-600 genes (same sets as iGEA),
BH FDR within each collection, FDR < 0.05 and overlap >= 3; each term linked to all its input genes.
Also 'ORA, all set sizes' (no 10-600 filter) at 10 collections. iGEA: terms kept (q < 0.10) in the 10-collection
CLI run, restricted to the collections of each step. Communities: network_modules.analyze_network_modules (as CLI).
Writes scaling_all.tsv and, per dataset, ora_terms_<dataset>_<set>.tsv / network_modules files for figures."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import glob, os, sys, logging
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection, exact_pvalue
from network_modules import analyze_network_modules, write_network_modules
H = os.path.dirname(os.path.abspath(__file__)); IN = os.path.join(H, "..", "inputs")
F = {"H: Hallmark Gene Sets": "h.all", "C2: Reactome Pathways": "c2.cp.reactome", "C2: KEGG Legacy": "c2.cp.kegg_legacy",
     "C2: KEGG MEDICUS": "c2.cp.kegg_medicus", "C2: BioCarta": "c2.cp.biocarta", "C2: Pathway Interaction Database": "c2.cp.pid",
     "C2: WikiPathways": "c2.cp.wikipathways", "C5: Gene Ontology: Biological Process": "c5.go.bp",
     "C5: Gene Ontology: Cellular Component": "c5.go.cc", "C5: Gene Ontology: Molecular Function": "c5.go.mf"}
GO = list(F)[-3:]
STEPS = [("GO (3)", GO), ("+ Reactome (4)", ["C2: Reactome Pathways"]), ("+ KEGG Legacy, KEGG MEDICUS (6)", ["C2: KEGG Legacy", "C2: KEGG MEDICUS"]),
         ("+ WikiPathways, BioCarta, PID (9)", ["C2: WikiPathways", "C2: BioCarta", "C2: Pathway Interaction Database"]),
         ("+ Hallmark (10)", ["H: Hallmark Gene Sets"])]
DATA = {"HIV": ("HIV_dependency_factors.symbols.txt", "HIV_dependency_factors.symbols_*"),
        "Hypoxia": ("hypoxia-genes.symbols.txt", "hypoxia-genes.symbols_*"),
        "Vaccination": ("vaccination_NAKAYA_HAI_28DY_3DY_POSITIVE.txt", "vaccination_NAKAYA_HAI_28DY_3DY_POSITIVE_*")}
bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter(); B = sorted(bg.genes)
libs = {n: GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n) for n, f in F.items()}
split = lambda s: [g.strip() for g in str(s).split(",") if g.strip()]

def ora(genes, lo, hi):
    out = {}
    for n in F:
        col = FastCollection(libs[n].library, B, n, lo, hi)
        x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in genes if g in col.index]; k = col.A @ x
        p = np.asarray(col._pvalues("Fisher's Exact Test", k, col.K, x.sum()), dtype=float)   # one-sided Fisher (vectorised)
        q = multipletests(p, method="fdr_bh")[1]
        out[n] = [dict(library=n, term=col.term_names[t], p=p[t], q=q[t],
                       genes=sorted(col.genes[j] for j in col.members[t] if col.genes[j] in genes))
                  for t in np.where((q < 0.05) & (k >= 3))[0]]
    return out

def stats(terms):
    deg = {}
    for t in terms:
        for g in t["genes"]: deg[g] = deg.get(g, 0) + 1
    m = analyze_network_modules(terms)
    return dict(terms=len(terms), genes=len(deg), edges=sum(deg.values()), max_gene_degree=max(deg.values()),
                communities=m.summary.get("n_communities"), modularity=m.summary.get("modularity")), m

rows = []
for ds, (fn, pat) in DATA.items():
    lst = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(os.path.join(IN, fn)) if l.strip()) if m]
    genes = set(GeneSet(lst, bg.genes, ds, hgcn=True, format=False).genes)
    o = ora(genes, 10, 600); o_all = ora(genes, 1, 10 ** 6)
    ig = pd.read_csv(glob.glob(f"{H}/{pat}/combined_iterative_results.tsv")[0], sep="\t"); ig = ig[ig["In network"] == True]
    igt = {n: [dict(library=n, term=r.Term, genes=split(r.Genes)) for r in ig[ig.Library == n].itertuples()] for n in F}
    cur = []
    for name, add in STEPS:
        cur += add
        for meth, src in (("ORA", o), ("iGEA", igt)):
            terms = [t for l in cur for t in src[l]]
            s, m = stats(terms)
            rows.append(dict(dataset=ds, input_genes=len(genes), collections=name, method=meth, **s))
            if name in ("GO (3)", "+ Hallmark (10)") and meth == "ORA":
                tag = "go" if name == "GO (3)" else "all"
                d = f"{H}/scaling_networks/{ds}_ORA_{tag}"; os.makedirs(d, exist_ok=True)
                pd.DataFrame([dict(Collection=t["library"], Term=t["term"], p=t["p"], FDR=t["q"], Genes=", ".join(t["genes"])) for t in terms]).to_csv(f"{d}/ora_terms.tsv", sep="\t", index=False)
                write_network_modules(m, d)
            if name in ("GO (3)", "+ Hallmark (10)") and meth == "iGEA":
                tag = "go" if name == "GO (3)" else "all"
                d = f"{H}/scaling_networks/{ds}_iGEA_{tag}"; os.makedirs(d, exist_ok=True)
                pd.DataFrame([dict(Collection=t["library"], Term=t["term"], Genes=", ".join(t["genes"])) for t in terms]).to_csv(f"{d}/igea_terms.tsv", sep="\t", index=False)
                write_network_modules(m, d)
    terms = [t for l in F for t in o_all[l]]; s, m = stats(terms)
    rows.append(dict(dataset=ds, input_genes=len(genes), collections="All 10, all set sizes", method="ORA (all set sizes)", **s))
    print(ds, "done", flush=True)
out = pd.DataFrame(rows); out.to_csv(f"{H}/scaling_all.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); print(out.to_string(index=False))
