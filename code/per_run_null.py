"""
Per-run random-list null distribution for iGEA networks.

For each analysis, ``n_lists`` random gene lists are drawn from the analysis' own background, with
exactly as many genes as the validated input list, and analysed with the same collections and the
same parameters (test, p-value threshold, minimum overlap, gene-set size limits, maximum number of
iterations). For every random list the cross-collection network is assembled and its largest
connected cluster is measured. The observed cluster is then compared with these values.

Statistics reported per metric:
  * z = (observed - mean) / SD of the random values;
  * empirical p = (r + 1) / (n + 1), where r is the number of random lists whose largest cluster
    is at least as large as the observed one (Phipson & Smyth 2010). With n = 1000 the smallest
    attainable value is 1/1001.

Random lists are reproducible: list i is drawn with numpy.random.default_rng([seed, size, i]) from
the alphabetically sorted background.
"""
import logging
import os
import time
from multiprocessing import get_context
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from fast_igea import FISHER, FastCollection
from network_connectivity_benchmark import NetworkConnectivityAnalyzer

logger = logging.getLogger(__name__)

DEFAULT_N_LISTS = 10000
DEFAULT_SEED = 2026

METRICS = [
    ("cluster_size", "Cluster Size"),
    ("n_genes", "Number of Genes"),
    ("n_terms", "Number of Terms"),
    ("n_edges", "Number of Edges"),
    ("avg_edges_per_gene", "Average Edges per Gene"),
]

_WORKER: Dict = {}


def largest_cluster_metrics(records: List[dict]) -> Dict[str, float]:
    """Largest connected cluster of the gene-term network built from iGEA records."""
    analyzer = NetworkConnectivityAnalyzer()
    analyzer.add_igea_results(records)
    clusters = analyzer.get_clusters()
    if not clusters:
        return {"cluster_size": 0, "n_genes": 0, "n_terms": 0, "n_edges": 0,
                "avg_edges_per_gene": 0.0, "n_libraries": 0, "n_terms_total": len(records)}
    c = clusters[0]
    return {"cluster_size": c["size"], "n_genes": c["n_genes"], "n_terms": c["n_terms"],
            "n_edges": c["n_edges"], "avg_edges_per_gene": c["avg_edges_per_gene"],
            "n_libraries": c["n_libraries"], "n_terms_total": len(records)}


def _init_worker(collections, background, size, seed, params):
    _WORKER.update(collections=collections, background=background, size=size, seed=seed, params=params)


def _one_list(i: int) -> Dict[str, float]:
    w = _WORKER
    rng = np.random.default_rng([w["seed"], w["size"], i])
    genes = [w["background"][j] for j in rng.choice(len(w["background"]), w["size"], replace=False)]
    records = []
    for col in w["collections"]:
        for it, term, p, removed in col.iterate(genes, **w["params"]):
            records.append({"Term": term, "Iteration": it, "Library": col.name, "iteration p-value": p,
                            "Genes removed for next iteration": removed})
    out = largest_cluster_metrics(records)
    out["list"] = i
    return out


def generate_null(libraries, background_genes, n_genes: int, p_value_method: str = FISHER,
                  p_threshold: float = 0.01, min_overlap: int = 3, min_term_size: int = 10,
                  max_term_size: int = 600, max_iterations: Optional[int] = None,
                  n_lists: int = DEFAULT_N_LISTS, seed: int = DEFAULT_SEED,
                  n_jobs: Optional[int] = None, progress_callback=None) -> pd.DataFrame:
    """Largest-cluster metrics for ``n_lists`` random lists matched to the analysis.

    Args:
        libraries: GeneSetLibrary objects used in the analysis.
        background_genes: the analysis background (the random lists are drawn from it).
        n_genes: size of the validated input list.
        n_jobs: worker processes (default: all CPUs; 1 = no multiprocessing).
    """
    background = sorted(set(background_genes))
    if n_genes > len(background):
        raise ValueError(f"Input size {n_genes} exceeds background size {len(background)}")
    collections = [FastCollection.from_library(lib, background, min_term_size, max_term_size) for lib in libraries]
    params = dict(p_threshold=p_threshold, min_overlap=min_overlap, max_iterations=max_iterations,
                  method=p_value_method)
    n_jobs = n_jobs or os.cpu_count() or 1
    t0 = time.time()
    if n_jobs == 1:
        _init_worker(collections, background, n_genes, seed, params)
        rows = []
        for i in range(n_lists):
            rows.append(_one_list(i))
            if progress_callback and (i + 1) % max(1, n_lists // 20) == 0:
                progress_callback(i + 1, n_lists)
    else:
        ctx = get_context("spawn")
        with ctx.Pool(n_jobs, initializer=_init_worker,
                      initargs=(collections, background, n_genes, seed, params)) as pool:
            rows = []
            for r in pool.imap(_one_list, range(n_lists), chunksize=max(1, n_lists // (n_jobs * 20))):
                rows.append(r)
                if progress_callback and len(rows) % max(1, n_lists // 20) == 0:
                    progress_callback(len(rows), n_lists)
    df = pd.DataFrame(rows).sort_values("list").reset_index(drop=True)
    df.attrs.update(n_genes=n_genes, n_lists=n_lists, seed=seed, n_jobs=n_jobs,
                    seconds=round(time.time() - t0, 1), background_size=len(background))
    logger.info(f"Per-run null: {n_lists} random lists of {n_genes} genes in {df.attrs['seconds']} s ({n_jobs} processes)")
    return df


def compare_to_null(observed: Dict[str, float], null: pd.DataFrame) -> pd.DataFrame:
    """z-score and empirical p-value for each largest-cluster metric."""
    rows = []
    n = len(null)
    for key, label in METRICS:
        x = null[key].to_numpy(dtype=float)
        obs = float(observed[key])
        sd = x.std()
        r = int((x >= obs).sum())
        rows.append({
            "Metric": label, "Observed": obs, "Null Mean": x.mean(), "Null SD": sd, "Null Max": x.max(),
            "Z-score": (obs - x.mean()) / sd if sd > 0 else float("nan"),
            "Random >= Observed": r, "Random Lists": n, "Empirical p": (r + 1) / (n + 1),
        })
    return pd.DataFrame(rows)


def format_report(table: pd.DataFrame, attrs: Dict, libraries: List[str], params: Dict) -> str:
    """Plain-text statistical report."""
    n = int(table["Random Lists"].iloc[0])
    lines = [
        "=" * 100,
        "iGEA NETWORK STATISTICAL REPORT (per-run random-list null)",
        "=" * 100,
        f"Random lists: {n}, each of {attrs['n_genes']} genes (= validated input size), drawn without replacement",
        f"from the analysis background ({attrs['background_size']} genes); seed {attrs['seed']}; {attrs['seconds']} s on {attrs['n_jobs']} process(es).",
        "Each random list was analysed with the same collections and parameters as the input:",
        "  " + ", ".join(f"{k}={v}" for k, v in params.items()),
        "Collections: " + "; ".join(libraries),
        "",
        "Statistics for the largest connected cluster:",
        "  z = (observed - random mean) / random SD",
        f"  empirical p = (r + 1) / (n + 1), r = random lists with a value >= observed, n = {n}",
        f"  (smallest attainable p = {1/(n+1):.4g}). The metrics describe the same cluster and are not independent tests.",
        "",
        f"  {'Metric':<26}{'Observed':>10}{'Rand mean':>11}{'Rand SD':>9}{'Rand max':>10}{'z':>8}{'r':>6}{'Emp. p':>10}",
        "  " + "-" * 90,
    ]
    for _, r in table.iterrows():
        lines.append(f"  {r['Metric']:<26}{r['Observed']:>10.4g}{r['Null Mean']:>11.3f}{r['Null SD']:>9.3f}"
                     f"{r['Null Max']:>10.4g}{r['Z-score']:>8.2f}{r['Random >= Observed']:>6d}{r['Empirical p']:>10.4g}")
    return "\n".join(lines) + "\n"
