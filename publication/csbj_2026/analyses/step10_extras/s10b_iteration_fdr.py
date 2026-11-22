"""(2) Chance selections by iteration index and by list size (R2-1, R5-M1).
Uses the cached 10,000 random lists of each dataset size (same lists as the paper statistics)."""
import os, numpy as np, pandas as pd
from common import *

def main():
    rows, size_rows = [], []
    for ds in ["HIV", "Hypoxia", "Vaccination"]:
        g = genes_of(ds); obs = records(g)
        cols = [FastCollection.from_library(l, BACKGROUND, 10, 600) for l in LIBS]
        params = dict(p_threshold=0.01, min_overlap=3, max_iterations=30, method=FISHER)
        gen = prn._generate_cached([(len(g), i) for i in range(10000)], cols, BACKGROUND, prn.DEFAULT_SEED, params, 4, None,
                                   os.path.expanduser("~/igea_cache_1s"))
        maxit = max([r[2] for r in obs] + [r[2] for L in gen for r in L])
        for k in range(1, maxit + 1):
            o = sum(1 for r in obs if r[2] == k); e = np.mean([sum(1 for r in L if r[2] == k) for L in gen])
            rows.append(dict(dataset=ds, input_genes=len(g), iteration=k, observed_terms=o, expected_chance_terms=round(e, 3),
                             chance_fraction=round(min(1, e / o), 3) if o else None))
        print(ds, "done", flush=True)
    pd.DataFrame(rows).to_csv("s10b_chance_by_iteration.tsv", sep="\t", index=False)
    # chance terms per random list as a function of list size (1,000 lists per size, no cache)
    cols = [FastCollection.from_library(l, BACKGROUND, 10, 600) for l in LIBS]
    prn._init_gen(cols, BACKGROUND, 2026, dict(p_threshold=0.01, min_overlap=3, max_iterations=30, method=FISHER))
    for n in [25, 50, 100, 200, 300, 400, 600, 800]:
        L = [prn._gen((n, i)) for i in range(1000)]
        cnt = [len(x) for x in L]; first = [sum(1 for r in x if r[2] == 1) for x in L]
        size_rows.append(dict(list_size=n, mean_chance_terms=round(np.mean(cnt), 2), p95=int(np.percentile(cnt, 95)),
                              mean_chance_terms_iteration1=round(np.mean(first), 2), max_iteration=max([r[2] for x in L for r in x] or [0])))
        print(size_rows[-1], flush=True)
    pd.DataFrame(size_rows).to_csv("s10b_chance_by_list_size.tsv", sep="\t", index=False)

if __name__ == '__main__':
    main()
