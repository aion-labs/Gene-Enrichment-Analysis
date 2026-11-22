"""Chunked version of per_run_null.run_statistics for runs longer than one shell call.
Identical random lists and computations (same _gen / ChanceReference / _score functions).
  gen   <gene_list> <work_dir> <task_start> <task_end> [n_lists]   -> work_dir/gen_<start>.pkl
  final <gene_list> <work_dir> <out_prefix> [n_lists] [q]           -> reports
Task order is the same as run_statistics: main lists (n_genes, 0..n_lists-1), then grid sizes."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, pickle, glob, logging, time
from multiprocessing import get_context
import numpy as np, pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection
import per_run_null as prn
P = [("H: Hallmark Gene Sets","h.all"),("C2: Reactome Pathways","c2.cp.reactome"),("C2: KEGG Legacy","c2.cp.kegg_legacy"),("C2: KEGG MEDICUS","c2.cp.kegg_medicus"),
     ("C2: BioCarta","c2.cp.biocarta"),("C2: Pathway Interaction Database","c2.cp.pid"),("C2: WikiPathways","c2.cp.wikipathways"),
     ("C5: Gene Ontology: Biological Process","c5.go.bp"),("C5: Gene Ontology: Cellular Component","c5.go.cc"),("C5: Gene Ontology: Molecular Function","c5.go.mf")]
PARAMS = dict(p_threshold=0.01, min_overlap=3, max_iterations=None, method="Fisher's Exact Test")
def setup(gl):
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    genes = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(gl) if l.strip()) if m]
    gs = GeneSet(genes, bg.genes, "x", hgcn=False, format=False)
    libs = [GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n) for n, f in P]
    return bg, gs, libs
def tasks_for(n, n_lists, q=0.10):
    grid = prn.grid_sizes(n); gl = max(200, n_lists // 10)
    return [(n, i) for i in range(n_lists)] + [(s, i) for s in grid if s != n for i in range(gl)], grid, gl
if __name__ == "__main__":
    mode, gl_path, work = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(work, exist_ok=True)
    bg, gs, libs = setup(gl_path); B = sorted(set(bg.genes))
    if mode == "gen":
        a, b = int(sys.argv[4]), int(sys.argv[5]); n_lists = int(sys.argv[6]) if len(sys.argv) > 6 else 10000
        tasks, _, _ = tasks_for(gs.size, n_lists); b = min(b, len(tasks))
        cols = [FastCollection.from_library(L, B) for L in libs]
        with get_context("spawn").Pool(os.cpu_count(), initializer=prn._init_gen, initargs=(cols, B, prn.DEFAULT_SEED, PARAMS)) as p:
            res = p.map(prn._gen, tasks[a:b], chunksize=50)
        pickle.dump((a, b, res), open(f"{work}/gen_{a:07d}.pkl", "wb")); print(f"tasks {a}-{b} of {len(tasks)} done")
    else:
        out = sys.argv[4]; n_lists = int(sys.argv[5]) if len(sys.argv) > 5 else 10000; qt = float(sys.argv[6]) if len(sys.argv) > 6 else 0.10
        tasks, grid, gl = tasks_for(gs.size, n_lists)
        gen = [None] * len(tasks)
        for f in sorted(glob.glob(f"{work}/gen_*.pkl")):
            a, b, res = pickle.load(open(f, "rb")); gen[a:b] = res
        assert all(g is not None for g in gen), "missing chunks"
        main = gen[:n_lists]; ref = prn.ChanceReference(len(libs)); ref.add(gs.size, main); pos = n_lists
        for s in grid:
            if s != gs.size: ref.add(s, gen[pos:pos + gl]); pos += gl
        names = [L.name for L in libs]
        with get_context("spawn").Pool(os.cpu_count(), initializer=prn._init_score, initargs=(names, ref, gs.size, qt)) as p:
            null = pd.DataFrame(p.map(prn._score, list(enumerate(main)), chunksize=50)).sort_values("list").reset_index(drop=True)
        obs = [(c, t, it, pp, g) for c, L in enumerate(libs) for it, t, pp, g in FastCollection.from_library(L, bg.genes).iterate(gs.genes, 0.01, 3, None)]
        q, m0 = prn.signal_adjusted_q(obs, gs.size, ref, qt); kept = [r for r, qq in zip(obs, q) if qq < qt]
        observed = prn.largest_cluster_metrics(prn._to_dicts(kept, names)); observed.update(n_terms_total=len(obs), n_terms_kept=len(kept))
        terms = pd.DataFrame([{"Library": names[c], "Iteration": it, "Term": t, "p-value": pp, "q-value": qq, "Kept": qq < qt, "Genes": ", ".join(g)} for (c, t, it, pp, g), qq in zip(obs, q)])
        table = prn.compare_to_null(observed, null)
        attrs = dict(n_genes=gs.size, n_lists=n_lists, seed=prn.DEFAULT_SEED, n_jobs=os.cpu_count(), background_size=len(B), q_threshold=qt, grid_sizes=grid,
                     grid_lists=gl, noise_genes=m0, expected_chance_terms=round(float(np.mean([len(m) for m in main])), 2), seconds="chunked")
        params = dict(p_value_method="Fisher's Exact Test", p_threshold=0.01, min_overlap=3, min_term_size=10, max_term_size=600, max_iterations=None, background="all_genes")
        open(out + "_report.txt", "w").write(prn.format_report(table, attrs, names, params, terms))
        terms.to_csv(out + "_terms.tsv", sep="\t", index=False); null.to_csv(out + "_random_lists.tsv", sep="\t", index=False)
        print({k: attrs[k] for k in ["noise_genes", "expected_chance_terms"]}, "kept", len(kept), "of", len(obs), {k: observed[k] for k in ["n_genes", "n_terms", "n_edges", "n_libraries"]})
        print(table[["Metric", "Observed", "Null Mean", "Null SD", "Z-score", "Random >= Observed", "Empirical p"]].round(4).to_string(index=False))
