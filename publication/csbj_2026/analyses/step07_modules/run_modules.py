"""Communities, hubs and bridges for the q<0.10 networks (HIV, hypoxia, vaccination, HIV GO-only).
Uses code/network_modules.py from the iterative_enrichment repo."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, pandas as pd
REPO = (str(_REPO) + "/code"); sys.path.insert(0, REPO)
from network_modules import analyze_network_modules, write_network_modules
S5 = "../step05_two_layer"
def load(f, go_only=False):
    d = pd.read_csv(f, sep="\t"); d = d[d["Kept"]]
    if go_only: d = d[d["Library"].str.startswith("C5")]
    return [dict(library=r.Library, term=r.Term, genes=[g.strip() for g in r.Genes.split(",")]) for r in d.itertuples()]
import glob
E9 = "../step09_end_to_end"
sets = {"hiv": load(glob.glob(f"{E9}/HIV_dependency_factors.symbols_*/HIV_dependency_factors.symbols_term_qvalues.tsv")[0]),
        "hiv_go": load("../../02 Figures/networks/hiv_go_terms.tsv"),
        "hypoxia": load(glob.glob(f"{E9}/hypoxia-genes.symbols_*/hypoxia-genes.symbols_term_qvalues.tsv")[0]),
        "vaccination": load(glob.glob(f"{E9}/vaccination_NAKAYA_HAI_28DY_3DY_POSITIVE_*/*_term_qvalues.tsv")[0])}
for name, t in sets.items():
    m = analyze_network_modules(t)
    write_network_modules(m, ".", prefix=name)
    print("=====", name, m.summary)
    print(m.communities[["Community","Terms","Genes","Collections","Lead term","Top hub genes (degree)","Bridge genes"]].to_string(index=False, max_colwidth=60))
