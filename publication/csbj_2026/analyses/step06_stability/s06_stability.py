"""Stability of iGEA under input perturbation (R1-1, R5-M5): 100 random 90 % subsamples of the HIV list.
All subsamples have the same size (238 genes = 90 % of 265), so they share one reference set of random lists (same background,
10 collections, same parameters): it is used for the q-values (two-level statistics, q < 0.10) and the network p.
Also: near ties in the full run (top two candidates within a factor 2 of p).
Usage: python s06_stability.py <out_dir> [n_ref_lists] [n_subsamples]"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging, itertools
from multiprocessing import get_context
import numpy as np, pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection, exact_pvalue
import per_run_null as prn
ALL = [("H: Hallmark Gene Sets","h.all"),("C2: Reactome Pathways","c2.cp.reactome"),("C2: KEGG Legacy","c2.cp.kegg_legacy"),("C2: KEGG MEDICUS","c2.cp.kegg_medicus"),
       ("C2: BioCarta","c2.cp.biocarta"),("C2: Pathway Interaction Database","c2.cp.pid"),("C2: WikiPathways","c2.cp.wikipathways"),
       ("C5: Gene Ontology: Biological Process","c5.go.bp"),("C5: Gene Ontology: Cellular Component","c5.go.cc"),("C5: Gene Ontology: Molecular Function","c5.go.mf")]
PARAMS = dict(p_threshold=0.01, min_overlap=3, max_iterations=30, method="Fisher's Exact Test")
W = {}
def init_sub(cols, names, ref, n):
    W.update(cols=cols, names=names, ref=ref, n=n)
def run_list(genes):
    recs = [(c, t, it, p, g) for c, col in enumerate(W["cols"]) for it, t, p, g in col.iterate(genes, **PARAMS)]
    q, _ = prn.signal_adjusted_q(recs, W["n"], W["ref"], 0.10)
    kept = [r for r, qq in zip(recs, q) if qq < 0.10]
    m = prn.largest_cluster_metrics(prn._to_dicts(kept, W["names"]))
    a = prn.NetworkConnectivityAnalyzer(); a.add_igea_results(prn._to_dicts(kept, W["names"])); cl = a.get_clusters()
    return {"selected": {r[1] for r in recs}, "kept": {r[1] for r in kept}, "cluster_genes": set(cl[0]["genes"]) if cl else set(),
            "cluster_terms": set(cl[0]["terms"]) if cl else set(), "metrics": m}
def jac(a, b): return len(a & b) / len(a | b) if (a | b) else 1.0
def near_ties(col, genes):
    """Iterations where the second-best candidate's p is within a factor 2 of the selected term's p."""
    x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in genes if g in col.index]; out = []; it = 1
    while True:
        n = x.sum(); k = col.A @ x; cand = np.flatnonzero(k >= 3)
        if n == 0 or cand.size == 0: break
        p = np.array([exact_pvalue("Fisher's Exact Test", int(k[c]), int(col.K[c]), int(n), col.M) for c in cand])
        o = np.argsort(p, kind="stable")
        if p[o[0]] >= 0.01: break
        second = p[o[1]] if len(o) > 1 else 1.0
        t = cand[o[0]]; out.append((it, col.term_names[t], p[o[0]], col.term_names[cand[o[1]]] if len(o) > 1 else "", second, second <= 2 * p[o[0]]))
        rem = col.members[t][x[col.members[t]] == 1]; x[rem] = 0; it += 1
    return out
if __name__ == "__main__":
    OUT = sys.argv[1]; NREF = int(sys.argv[2]) if len(sys.argv) > 2 else 5000; NSUB = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter(); B = sorted(bg.genes)
    hiv = sorted(GeneSet([m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m], bg.genes, "HIV", hgcn=True, format=False).genes)  # one universe
    libs = [GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n) for n, f in ALL]; names = [n for n, _ in ALL]
    cols = [FastCollection.from_library(L, B) for L in libs]
    n = int(round(0.9 * len(hiv)))
    grid = prn.grid_sizes(n); GL = max(200, NREF // 10)
    tasks = [(n, i) for i in range(NREF)] + [(s, i) for s in grid if s != n for i in range(GL)]
    ctx = get_context("spawn")
    gen = prn._generate_cached(tasks, cols, B, prn.DEFAULT_SEED, PARAMS, os.cpu_count(), None,
                               os.path.expanduser("~/igea_cache"))   # resumable
    ref = prn.ChanceReference(len(cols)); ref.add(n, gen[:NREF]); pos = NREF
    for s in grid:
        if s != n: ref.add(s, gen[pos:pos + GL]); pos += GL
    with ctx.Pool(os.cpu_count(), initializer=prn._init_score, initargs=(names, ref, n, 0.10)) as p:
        null = pd.DataFrame(p.map(prn._score, list(enumerate(gen[:NREF])), chunksize=50))
    rng = np.random.default_rng(2026)
    subs = [sorted(rng.choice(hiv, n, replace=False)) for _ in range(NSUB)]
    # full list scored against its own reference (as in the paper)
    import glob
    full = pd.read_csv(glob.glob(os.path.join(os.path.dirname(OUT), "step09_end_to_end/HIV_dependency_factors.symbols_*/HIV_dependency_factors.symbols_term_qvalues.tsv"))[0], sep="\t")
    full_kept = set(full[full.Kept.astype(bool)].Term); full_sel = set(full.Term)
    with ctx.Pool(os.cpu_count(), initializer=init_sub, initargs=(cols, names, ref, n)) as p:
        res = p.map(run_list, subs, chunksize=5)
    # full-list cluster genes from the paper network
    fg = set(); [fg.update(str(g).split(", ")) for g in full[full.Kept.astype(bool)].Genes]
    rows = []
    for i, r in enumerate(res):
        m = r["metrics"]; x = null["n_genes"].to_numpy()
        rows.append({"subsample": i, "selected": len(r["selected"]), "kept": len(r["kept"]), "jaccard_selected_vs_full": jac(r["selected"], full_sel),
                     "jaccard_kept_vs_full": jac(r["kept"], full_kept), "full_kept_recovered": len(r["kept"] & full_kept), "new_kept_terms": len(r["kept"] - full_kept),
                     "cluster_genes": m["n_genes"], "cluster_terms": m["n_terms"], "cluster_gene_jaccard_vs_full": jac(r["cluster_genes"], fg),
                     "network_p_genes": ((x >= m["n_genes"]).sum() + 1) / (len(x) + 1)})
    d = pd.DataFrame(rows); d.to_csv(f"{OUT}/s06_subsamples.tsv", sep="\t", index=False)
    rec = pd.DataFrame([{"Library": r.Library, "Term": r.Term, "Iteration": r.Iteration, "q_full": r["q-value"],
                         "kept_full": bool(r.Kept), "recurrence_selected": np.mean([r.Term in x["selected"] for x in res]),
                         "recurrence_kept": np.mean([r.Term in x["kept"] for x in res])} for _, r in full.iterrows()])
    rec.to_csv(f"{OUT}/s06_term_recurrence.tsv", sep="\t", index=False)
    nt = [(name,) + t for (name, _), col in zip(ALL, cols) for t in near_ties(col, hiv)]
    nt = pd.DataFrame(nt, columns=["Library", "Iteration", "Selected", "p", "Runner-up", "p runner-up", "near tie (<=2x)"])
    nt = nt.merge(rec[["Term", "recurrence_kept", "kept_full"]], left_on="Selected", right_on="Term", how="left").drop(columns="Term")
    nt.to_csv(f"{OUT}/s06_near_ties.tsv", sep="\t", index=False)
    print(d.describe().loc[["mean", "std", "min", "50%", "max"]].round(3).T.to_string())
    k = rec[rec.kept_full]
    print(f"\nFull-list kept terms ({len(k)}): recurrence among kept in subsamples: median {k.recurrence_kept.median():.2f}; >=0.8: {(k.recurrence_kept>=0.8).sum()}; >=0.5: {(k.recurrence_kept>=0.5).sum()}; <0.5: {(k.recurrence_kept<0.5).sum()}")
    print(f"Near ties in full run: {int(nt['near tie (<=2x)'].sum())} of {len(nt)} iterations; recurrence (kept) of near-tie terms: median {nt[nt['near tie (<=2x)']].recurrence_kept.median():.2f} vs others {nt[~nt['near tie (<=2x)']].recurrence_kept.median():.2f}")
    print("Subsamples with network p <= 0.001 (genes):", int((d.network_p_genes <= 0.001).sum()), "of", len(d))
