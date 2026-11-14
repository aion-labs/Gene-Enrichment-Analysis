"""
Streamlit section: network statistics against random gene lists generated for this run.

The random lists have the validated input size, are drawn from the run's background and are
analysed with the run's collections and parameters (see per_run_null.py). Computation starts only
when the user presses the button, so the page stays responsive.
"""
import logging
import time
from typing import Callable, Dict, Optional

import pandas as pd

from per_run_null import (DEFAULT_SEED, compare_to_null, format_report, generate_null,
                          largest_cluster_metrics)

logger = logging.getLogger(__name__)


def _run_signature(iter_enrich: Dict) -> tuple:
    first = next(iter(iter_enrich.values()))
    return (tuple(sorted(first.gene_set.genes)), tuple(sorted(iter_enrich.keys())),
            first.background_gene_set.name, first.p_value_method_name, first.p_threshold,
            first.min_overlap, first.min_term_size, first.max_term_size, first.max_iterations,
            tuple((k, len(v.results)) for k, v in sorted(iter_enrich.items())))


def render_per_run_benchmark(st, state, download_link: Callable, default_n_lists: int = 10000,
                             max_n_lists: int = 100000, n_jobs: Optional[int] = None,
                             seconds_per_list_hint: Optional[float] = None) -> None:
    """Render the statistics section for the current iGEA results."""
    iter_enrich = state.get("iter_enrich") or {}
    if not iter_enrich:
        st.info("Run iterative enrichment analysis first to compute network statistics.")
        return

    first = next(iter(iter_enrich.values()))
    sig = _run_signature(iter_enrich)
    stored = state.get("per_run_benchmark")
    if stored and stored.get("signature") != sig:
        state.per_run_benchmark = stored = None

    n_genes = first.gene_set.size
    st.caption(
        f"Your network is compared with networks built from random gene lists of exactly {n_genes} genes "
        f"(your validated input), drawn from your background ({first.background_gene_set.size} genes) and "
        "analysed with the same collections and parameters."
    )
    c1, c2 = st.columns(2)
    n_lists = c1.number_input("Number of random lists", min_value=100, max_value=max_n_lists,
                              value=default_n_lists, step=100, key="null_n_lists")
    seed = c2.number_input("Random seed", min_value=0, value=DEFAULT_SEED, step=1, key="null_seed")
    if seconds_per_list_hint:
        st.caption(f"Estimated time: about {max(1, round(n_lists * seconds_per_list_hint / 60))} min "
                   "(longer for large gene lists).")

    if st.button("Compute network statistics", type="primary"):
        records = []
        for lib_name, it in iter_enrich.items():
            for r in it.results:
                records.append({"Term": r.get("Term_raw", r.get("Term")), "Iteration": r.get("Iteration"),
                                "Library": lib_name,
                                "Genes removed for next iteration": r.get("Genes removed for next iteration", [])})
        observed = largest_cluster_metrics(records)
        if observed["n_genes"] == 0:
            st.warning("No network was formed, so there is nothing to compare.")
            return
        params = dict(p_value_method=first.p_value_method_name, p_threshold=first.p_threshold,
                      min_overlap=first.min_overlap, min_term_size=first.min_term_size,
                      max_term_size=first.max_term_size, max_iterations=first.max_iterations)
        bar = st.progress(0.0, text="Analysing random gene lists…")
        t0 = time.time()

        def progress(done, total):
            bar.progress(done / total, text=f"Analysing random gene lists… {done:,}/{total:,} "
                                            f"({time.time() - t0:.0f} s)")
        libraries = [it.gene_set_library for it in iter_enrich.values()]
        null = generate_null(libraries, first.background_gene_set.genes, n_genes, n_lists=int(n_lists),
                             seed=int(seed), n_jobs=n_jobs, progress_callback=progress, **params)
        bar.empty()
        table = compare_to_null(observed, null)
        shown = dict(params, background=first.background_gene_set.name)
        name = state.get("gene_set_name") or "gene_list"
        report = (f"Gene list: {name}\nValidated input genes: {n_genes}\n"
                  f"Largest cluster: {observed['cluster_size']} nodes, {observed['n_genes']} genes, "
                  f"{observed['n_terms']} terms, {observed['n_edges']} edges, "
                  f"{observed['avg_edges_per_gene']:.4f} edges per gene, {observed['n_libraries']} collections\n\n"
                  + format_report(table, null.attrs, [lib.name for lib in libraries], shown))
        state.per_run_benchmark = stored = {"signature": sig, "table": table, "report": report,
                                            "null": null, "observed": observed, "attrs": dict(null.attrs)}

    if not stored:
        return
    obs, attrs, table = stored["observed"], stored["attrs"], stored["table"]
    st.success(f"Largest cluster: {obs['n_genes']} genes, {obs['n_terms']} terms, {obs['n_edges']} edges, "
               f"{obs['n_libraries']} collections. Compared with {attrs['n_lists']:,} random lists "
               f"(seed {attrs['seed']}, {attrs['seconds']} s).")
    show = table[["Metric", "Observed", "Null Mean", "Null SD", "Z-score", "Random >= Observed", "Empirical p"]].copy()
    show.columns = ["Metric", "Your network", "Random mean", "Random SD", "z-score", "Random lists ≥ yours", "Empirical p"]
    st.dataframe(show.round(4), use_container_width=True, hide_index=True)
    with st.expander("📖 How to read these numbers"):
        n = attrs["n_lists"]
        st.markdown(f"""
- **Empirical p** = (r + 1) / (n + 1), where *r* is the number of random lists whose largest cluster
  was at least as large as yours and *n* = {n:,}. The smallest possible value is {1/(n+1):.2g}.
  This makes no assumption about the shape of the random distribution.
- **z-score** = (your value − random mean) / random SD, shown for reference. The random values are
  counts with skewed distributions, so read the empirical p rather than converting z to a probability.
- The five metrics describe the same cluster and are not independent tests.
- A small p means your genes form a more connected cross-collection network than random lists of the
  same size, background and settings. It does not test individual terms.
""")
    name = state.get("gene_set_name") or "statistical_report"
    st.markdown(f"📄 {download_link(stored['report'], name + '_statistical_report', 'txt')} · "
                f"{download_link(stored['null'].to_csv(sep=chr(9), index=False), name + '_random_lists', 'tsv')}",
                unsafe_allow_html=True)
