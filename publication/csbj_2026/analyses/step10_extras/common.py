"""Shared setup for the step10 reviewer analyses (same data, background, collections and defaults as the CLI)."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import os, sys, logging
import numpy as np
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection, FISHER
import per_run_null as prn
H = os.path.dirname(os.path.abspath(__file__)); IN = os.path.join(H, "..", "inputs")
ORDER = [("H: Hallmark Gene Sets", "h.all"), ("C2: BioCarta", "c2.cp.biocarta"), ("C2: KEGG MEDICUS", "c2.cp.kegg_medicus"),
         ("C2: Pathway Interaction Database", "c2.cp.pid"), ("C2: Reactome Pathways", "c2.cp.reactome"),
         ("C2: WikiPathways", "c2.cp.wikipathways"), ("C2: KEGG Legacy", "c2.cp.kegg_legacy"),
         ("C5: Gene Ontology: Biological Process", "c5.go.bp"), ("C5: Gene Ontology: Cellular Component", "c5.go.cc"),
         ("C5: Gene Ontology: Molecular Function", "c5.go.mf")]
DATA = {"HIV": "HIV_dependency_factors.symbols.txt", "Hypoxia": "hypoxia-genes.symbols.txt",
        "Vaccination": "vaccination_NAKAYA_HAI_28DY_3DY_POSITIVE.txt",
        "Random control": os.path.join("..", "step09_end_to_end", "inputs", "random_control_seed999_265genes.txt")}
bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); BACKGROUND = sorted(set(bg.genes))
LIBS = [GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n) for n, f in ORDER]
DEFAULT = dict(p_threshold=0.01, min_overlap=3, min_term_size=10, max_term_size=600, max_iterations=30)
_conv = GeneConverter()

def genes_of(name):
    lst = [m for m in (_conv.validate_and_map_symbol(l.strip()) for l in open(os.path.join(IN, DATA[name])) if l.strip()) if m]
    return sorted(set(GeneSet(lst, bg.genes, name, hgcn=True, format=False).genes))

def records(genes, p_threshold=0.01, min_overlap=3, min_term_size=10, max_term_size=600, max_iterations=30):
    out = []
    for c, lib in enumerate(LIBS):
        col = FastCollection.from_library(lib, BACKGROUND, min_term_size, max_term_size)
        for it, term, p, removed in col.iterate(genes, p_threshold=p_threshold, min_overlap=min_overlap,
                                                 max_iterations=max_iterations, method=FISHER):
            out.append((c, term, it, p, removed))
    return out

def stats(genes, n_lists=10000, n_jobs=4, cache_dir=None, **params):
    p = {**DEFAULT, **params}
    recs = records(genes, **p)
    table, null, terms, observed, attrs = prn.run_statistics(recs, LIBS, BACKGROUND, len(genes), p_value_method=FISHER,
                                                            n_lists=n_lists, n_jobs=n_jobs, cache_dir=cache_dir, **p)
    return recs, table, null, terms, observed, attrs
