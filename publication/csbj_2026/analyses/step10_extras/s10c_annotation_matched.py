"""(3) Annotation-matched random lists (R5-M2, R2-1): is network significance explained by annotation bias?
Each random list draws, for every input gene, a background gene with a similar annotation count (number of the
10,439 gene sets of 10-600 genes that contain it; bin 0 plus nine quantile bins). Reduced-size lists for the q
adjustment are random subsets of matched lists. Scoring, q-values and network test exactly as run_statistics."""
import sys
from multiprocessing import get_context
import numpy as np, pandas as pd
from common import *
import matched_lib

def main(ds, n_lists=10000):
    cols = [FastCollection.from_library(l, BACKGROUND, 10, 600) for l in LIBS]
    ann = {g: 0 for g in BACKGROUND}
    for col in cols:
        for t in range(len(col.term_names)):
            for j in col.members[t]: ann[col.genes[j]] += 1
    vals = np.array([ann[g] for g in BACKGROUND]); pos = vals[vals > 0]
    edges = np.unique(np.quantile(pos, np.linspace(0, 1, 10)))
    bin_of = lambda a: 0 if a == 0 else 1 + int(np.searchsorted(edges[1:-1], a, side="right"))
    bins = {}
    for g in BACKGROUND: bins.setdefault(bin_of(ann[g]), []).append(g)
    genes = genes_of(ds); n = len(genes); need = {}
    for g in genes: need[bin_of(ann[g])] = need.get(bin_of(ann[g]), 0) + 1
    params = dict(p_threshold=0.01, min_overlap=3, max_iterations=30, method=FISHER)
    grid = prn.grid_sizes(n); gl = max(200, n_lists // 10)
    tasks = [(n, i) for i in range(n_lists)] + [(s, i) for s in grid if s != n for i in range(gl)]
    with get_context("spawn").Pool(4, initializer=matched_lib.init, initargs=(cols, BACKGROUND, need, bins, n, params)) as pool:
        gen = pool.map(matched_lib.gen, tasks, chunksize=200)
    main_l = gen[:n_lists]; ref = prn.ChanceReference(len(cols)); ref.add(n, main_l); p = n_lists
    for s in grid:
        if s != n: ref.add(s, gen[p:p + gl]); p += gl
    names = [l.name for l in LIBS]
    scored = prn._pool_map(prn._score, list(enumerate(main_l)), prn._init_score, (names, ref, n, 0.10), 4)
    null = pd.DataFrame(scored).sort_values("list").reset_index(drop=True)
    obs = records(genes)
    q, m0 = prn.signal_adjusted_q(obs, n, ref, 0.10)
    kept = [r for r, qq in zip(obs, q) if qq < 0.10]
    observed = prn.largest_cluster_metrics(prn._to_dicts(kept, names))
    table = prn.compare_to_null(observed, null)
    table.insert(0, "null", "annotation-matched"); table.insert(0, "dataset", ds)
    table["terms_kept"] = len(kept); table["terms_selected"] = len(obs)
    table["chance_terms_per_list"] = round(float(np.mean([len(m) for m in main_l])), 2)
    table["input_mean_annotation"] = round(float(np.mean([ann[g] for g in genes])), 1)
    table["matched_lists_mean_annotation"] = round(float(np.mean([np.mean([ann[g] for g in []]) if False else 0])), 1)
    table["background_mean_annotation"] = round(float(vals.mean()), 1)
    tag = ds.replace(" ", "_")
    table.to_csv(f"s10c_matched_{tag}.tsv", sep="\t", index=False)
    pd.DataFrame([dict(Library=names[c], Term=t, Iteration=it, p=pp, q=qq, Kept=qq < 0.10) for (c, t, it, pp, _), qq in zip(obs, q)]).to_csv(
        f"s10c_matched_{tag}_terms.tsv", sep="\t", index=False)
    print(ds, "kept", len(kept), "of", len(obs), "| chance/list", table.chance_terms_per_list[0],
          "| input annotation", table.input_mean_annotation[0], "vs background", table.background_mean_annotation[0])
    print(table.iloc[:, 2:9].to_string(index=False), flush=True)

if __name__ == "__main__":
    main(sys.argv[1])
