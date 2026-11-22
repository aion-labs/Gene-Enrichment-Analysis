"""(4) Uncertainty of the empirical statistics (R2-1): bootstrap intervals for term q-values (resampling the
random lists, 200 replicates) and exact binomial intervals for the network empirical p-values.
Uses the cached random lists of the paper analysis (identical tasks to run_statistics)."""
import os, numpy as np, pandas as pd
from scipy.stats import beta
from common import *

def main():
    rows, prow = [], []
    for ds in ["HIV", "Hypoxia", "Vaccination"]:
        g = genes_of(ds); n = len(g); obs = records(g)
        cols = [FastCollection.from_library(l, BACKGROUND, 10, 600) for l in LIBS]
        params = dict(p_threshold=0.01, min_overlap=3, max_iterations=30, method=FISHER)
        grid = prn.grid_sizes(n); gl = max(200, 10000 // 10)
        tasks = [(n, i) for i in range(10000)] + [(s, i) for s in grid if s != n for i in range(gl)]
        gen = prn._generate_cached(tasks, cols, BACKGROUND, prn.DEFAULT_SEED, params, 4, None, os.path.expanduser("~/igea_cache_1s"))
        main_l = gen[:10000]; blocks = {}; pos = 10000
        for s in grid:
            if s != n: blocks[s] = gen[pos:pos + gl]; pos += gl
        def q_of(idx_main, idx_grid):
            ref = prn.ChanceReference(len(cols)); ref.add(n, [main_l[i] for i in idx_main])
            for s, L in blocks.items(): ref.add(s, [L[i] for i in idx_grid[s]])
            return prn.signal_adjusted_q(obs, n, ref, 0.10)[0]
        q0 = q_of(range(10000), {s: range(gl) for s in blocks})
        rng = np.random.default_rng(2026); boots = []
        for b in range(200):
            boots.append(q_of(rng.integers(0, 10000, 10000), {s: rng.integers(0, gl, gl) for s in blocks}))
        B = np.array(boots)
        for k, r in enumerate(obs):
            lo, hi = np.percentile(B[:, k], [2.5, 97.5])
            rows.append(dict(dataset=ds, collection=LIBS[r[0]].name, term=r[1], iteration=r[2], p=r[3], q=round(q0[k], 4),
                             q_low=round(lo, 4), q_high=round(hi, 4), kept=q0[k] < 0.10,
                             kept_in_bootstrap_fraction=round(float((B[:, k] < 0.10).mean()), 3)))
        print(ds, "kept", int((q0 < 0.1).sum()), "of", len(obs), "; terms whose 95% interval crosses 0.10:",
              sum(1 for k in range(len(obs)) if np.percentile(B[:, k], 2.5) < 0.1 <= np.percentile(B[:, k], 97.5)), flush=True)
    pd.DataFrame(rows).to_csv("s10d_q_bootstrap.tsv", sep="\t", index=False)
    # exact (Clopper-Pearson) 95% intervals for r/n, network metrics of the paper reports
    for ds, r_list in [("HIV", {"Cluster Size": 0, "Average Edges per Gene": 62}), ("HIV GO (Figure 2)", {"Cluster Size": 0, "Average Edges per Gene": 46})]:
        for m, r in r_list.items():
            nn = 10000; lo = 0 if r == 0 else beta.ppf(0.025, r, nn - r + 1); hi = beta.ppf(0.975, r + 1, nn - r)
            prow.append(dict(network=ds, metric=m, r=r, n=nn, empirical_p=round((r + 1) / (nn + 1), 5), tail_prob_low=round(lo, 5), tail_prob_high=round(hi, 5)))
    pd.DataFrame(prow).to_csv("s10d_network_p_intervals.tsv", sep="\t", index=False)
    print(pd.DataFrame(prow).to_string(index=False))

if __name__ == "__main__":
    main()
