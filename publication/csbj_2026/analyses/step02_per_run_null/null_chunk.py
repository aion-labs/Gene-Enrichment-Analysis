"""Per-run null in chunks (identical to per_run_null.generate_null: list i uses
default_rng([seed, size, i]) on the sorted background), for runs longer than one shell call.
Usage: python null_chunk.py <gene_list> <start> <end> <out_tsv> [seed]"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os
from multiprocessing import get_context
import pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code")
import logging; logging.disable(logging.WARNING)
import per_run_null as prn
from fast_igea import FastCollection
from background_gene_set import BackgroundGeneSet
from gene_converter import GeneConverter
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
LIBS = {"H: Hallmark Gene Sets": "h.all", "C2: BioCarta": "c2.cp.biocarta", "C2: KEGG MEDICUS": "c2.cp.kegg_medicus",
        "C2: Pathway Interaction Database": "c2.cp.pid", "C2: Reactome Pathways": "c2.cp.reactome",
        "C2: WikiPathways": "c2.cp.wikipathways", "C2: KEGG Legacy": "c2.cp.kegg_legacy",
        "C5: Gene Ontology: Biological Process": "c5.go.bp", "C5: Gene Ontology: Cellular Component": "c5.go.cc",
        "C5: Gene Ontology: Molecular Function": "c5.go.mf"}
if __name__ == "__main__":
    gl, start, end, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else prn.DEFAULT_SEED
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt")          # validated, as in the CLI
    conv = GeneConverter()
    genes = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(gl) if l.strip()) if m]
    gs = GeneSet(genes, bg.genes, "x", hgcn=False, format=False)
    background = sorted(set(bg.genes))
    cols = [FastCollection.from_library(GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n), background) for n, f in LIBS.items()]
    params = dict(p_threshold=0.01, min_overlap=3, max_iterations=None, method="Fisher's Exact Test")
    with get_context("spawn").Pool(os.cpu_count(), initializer=prn._init_worker, initargs=(cols, background, gs.size, seed, params)) as p:
        rows = p.map(prn._one_list, range(start, end), chunksize=20)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False); print(f"lists {start}-{end} of size {gs.size} -> {out}")
