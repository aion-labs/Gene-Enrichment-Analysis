"""Empirical p-values for the HIV network, using exactly the null the iGEA code builds
(parallel_null_distribution.compute_null_distribution_from_raw_permutations: size bin 300,
10 collections, p <= 0.01, iteration <= 30). Per-run values are captured by wrapping the
DataFrame the function builds, so z-scores and empirical p use identical numbers.
Run: python s01c_empirical_p.py <repo> <out_dir>"""
import sys, json, logging
from pathlib import Path
import numpy as np, pandas as pd
REPO = Path(sys.argv[1]).resolve(); OUT = Path(sys.argv[2]).resolve()
sys.path.insert(0, str(REPO / "code")); logging.disable(logging.WARNING)
import parallel_null_distribution as pnd
captured = {}
_orig = pnd.pd.DataFrame
class _Capture(_orig):
    def __init__(self, data=None, *a, **k):
        super().__init__(data, *a, **k)
        if isinstance(data, list) and data and isinstance(data[0], dict) and "largest_cluster_genes" in data[0]:
            captured["runs"] = _orig(data)
pnd.pd.DataFrame = _Capture
LIBS = ["H: Hallmark Gene Sets","C2: CP: BioCarta","C2: CP: KEGG MEDICUS","C2: CP: Pathway Interaction Database",
        "C2: CP: Reactome Pathways","C2: CP: WikiPathways","C2: CP: KEGG Legacy",
        "C5: GO: Biological Process","C5: GO: Cellular Component","C5: GO: Molecular Function"]
stats, size = pnd.compute_null_distribution_from_raw_permutations(REPO / "permutations/merged_permutation_results.tsv", 272, LIBS, 0.01, 30)
pnd.pd.DataFrame = _orig
runs = captured["runs"]; runs.to_csv(OUT / "s01_null_runs_size300.csv", index=False)
obs = {"largest_component_size": ("Cluster size (genes + terms)", 241), "largest_cluster_genes": ("Genes", 179),
       "largest_cluster_terms": ("Terms", 62), "largest_cluster_edges": ("Edges", 412),
       "largest_cluster_avg_edges_per_gene": ("Edges per gene", 412 / 179)}
rows = []
for m, (label, v) in obs.items():
    x = runs[m].to_numpy(float); r = int((x >= v).sum()); n = len(x)
    rows.append({"Metric": label, "Observed": round(v, 4), "Null mean": round(x.mean(), 2), "Null SD": round(x.std(), 2),
                 "Null max": round(x.max(), 2), "z": round((v - x.mean()) / x.std(), 2),
                 "Random >= observed (r)": r, "n random": n, "Empirical p = (r+1)/(n+1)": round((r + 1) / (n + 1), 4),
                 "Null skewness": round(float(pd.Series(x).skew()), 2),
                 "check: code null mean": round(stats[m]["mean"], 4)})
t = pd.DataFrame(rows); t.to_csv(OUT / "s01_hiv_z_and_empirical_p.csv", index=False); print("null size bin:", size); print(t.to_string(index=False))
