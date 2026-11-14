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

from per_run_null import (DEFAULT_GRID_LISTS, DEFAULT_Q_THRESHOLD, DEFAULT_SEED, format_report,
                          run_statistics)

logger = logging.getLogger(__name__)


def _run_signature(iter_enrich: Dict) -> tuple:
    first = next(iter(iter_enrich.values()))
    return (tuple(sorted(first.gene_set.genes)), tuple(sorted(iter_enrich.keys())),
            first.background_gene_set.name, first.p_value_method_name, first.p_threshold,
            first.min_overlap, first.min_term_size, first.max_term_size, first.max_iterations,
            tuple((k, len(v.results)) for k, v in sorted(iter_enrich.items())))


def render_per_run_benchmark(st, state, download_link: Callable, default_n_lists: int = 10000,
                             max_n_lists: int = 100000, n_jobs: Optional[int] = None, grid_lists: Optional[int] = None,
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
    c1, c2, c3 = st.columns(3)
    n_lists = c1.number_input("Number of random lists", min_value=100, max_value=max_n_lists,
                              value=default_n_lists, step=100, key="null_n_lists")
    q_threshold = c2.number_input("Term q-value threshold (0 = keep all)", min_value=0.0, max_value=1.0,
                                  value=DEFAULT_Q_THRESHOLD, step=0.01, format="%.2f", key="null_q_threshold")
    seed = c3.number_input("Random seed", min_value=0, value=DEFAULT_SEED, step=1, key="null_seed")
    if seconds_per_list_hint:
        st.caption(f"Estimated time: about {max(1, round(n_lists * seconds_per_list_hint / 60))} min "
                   "(longer for large gene lists).")

    if st.button("Compute network statistics", type="primary"):
        libraries = [it.gene_set_library for it in iter_enrich.values()]
        records = []
        for c, it in enumerate(iter_enrich.values()):
            for r in it.results:
                records.append((c, r.get("Term_raw", r.get("Term")), r.get("Iteration"), r.get("iteration p-value"),
                                list(r.get("Genes removed for next iteration", []))))
        if not records:
            st.warning("No terms were selected, so there is nothing to compare.")
            return
        params = dict(p_value_method=first.p_value_method_name, p_threshold=first.p_threshold,
                      min_overlap=first.min_overlap, min_term_size=first.min_term_size,
                      max_term_size=first.max_term_size, max_iterations=first.max_iterations)
        bar = st.progress(0.0, text="Analysing random gene lists…")
        t0 = time.time()

        def progress(done, total, label):
            bar.progress(done / total, text=f"Analysing {label}… {done:,}/{total:,} ({time.time() - t0:.0f} s)")
        qt = float(q_threshold) if q_threshold and q_threshold > 0 else None
        table, null, terms, observed, attrs = run_statistics(
            records, libraries, first.background_gene_set.genes, n_genes, n_lists=int(n_lists), seed=int(seed),
            q_threshold=qt, grid_lists=grid_lists, n_jobs=n_jobs, progress_callback=progress, **params)
        bar.empty()
        shown = dict(params, background=first.background_gene_set.name)
        name = state.get("gene_set_name") or "gene_list"
        report = (f"Gene list: {name}\nValidated input genes: {n_genes}\n"
                  f"Largest cluster: {observed['cluster_size']} nodes, {observed['n_genes']} genes, "
                  f"{observed['n_terms']} terms, {observed['n_edges']} edges, "
                  f"{observed['avg_edges_per_gene']:.4f} edges per gene, {observed['n_libraries']} collections\n\n"
                  + format_report(table, attrs, [lib.name for lib in libraries], shown, terms))
        state.per_run_benchmark = stored = {"signature": sig, "table": table, "report": report, "terms": terms,
                                            "null": null, "observed": observed, "attrs": dict(attrs)}

    if not stored:
        return
    obs, attrs, table = stored["observed"], stored["attrs"], stored["table"]
    qt = attrs.get("q_threshold")
    if qt is not None:
        st.info(f"**Term level:** {obs['n_terms_kept']} of {obs['n_terms_total']} terms have q < {qt} and are kept. "
                f"A random list of this size selects about {attrs['expected_chance_terms']} terms by chance; "
                f"random lists are filtered the same way before their networks are built.")
    st.success(f"Largest cluster: {obs['n_genes']} genes, {obs['n_terms']} terms, {obs['n_edges']} edges, "
               f"{obs['n_libraries']} collections. Compared with {attrs['n_lists']:,} random lists "
               f"(seed {attrs['seed']}, {attrs['seconds']} s).")
    show = table[["Metric", "Observed", "Null Mean", "Null SD", "Z-score", "Random >= Observed", "Empirical p"]].copy()
    show.columns = ["Metric", "Your network", "Random mean", "Random SD", "z-score", "Random lists ≥ yours", "Empirical p"]
    st.dataframe(show.round(4), use_container_width=True, hide_index=True)
    if "terms" in stored and len(stored["terms"]):
        with st.expander("Term q-values"):
            st.dataframe(stored["terms"][["Library", "Iteration", "Term", "p-value", "q-value", "Kept"]],
                         use_container_width=True, hide_index=True)
    with st.expander("📖 How to read these numbers"):
        n = attrs["n_lists"]
        st.markdown(f"""
- **Empirical p** = (r + 1) / (n + 1), where *r* is the number of random lists whose largest cluster
  was at least as large as yours and *n* = {n:,}. The smallest possible value is {1/(n+1):.2g}.
  This makes no assumption about the shape of the random distribution.
- **z-score** = (your value − random mean) / random SD, shown for reference. The random values are
  counts with skewed distributions, so read the empirical p rather than converting z to a probability.
- The five metrics describe the same cluster and are not independent tests.
- **q-value (per term)**: estimated fraction of chance terms among terms at least this strong, from
  random lists of the same size, background and settings, adjusted for genes already explained by
  accepted terms. Terms above the threshold are removed from your network and from every random network.
- A small network p means your genes form a more connected cross-collection network than random lists
  of the same size, background and settings.
""")
    name = state.get("gene_set_name") or "statistical_report"
    st.markdown(f"📄 {download_link(stored['report'], name + '_statistical_report', 'txt')} · "
                f"{download_link(stored['null'].to_csv(sep=chr(9), index=False), name + '_random_lists', 'tsv')}"
                + (f" · {download_link(stored['terms'].to_csv(sep=chr(9), index=False), name + '_term_qvalues', 'tsv')}"
                   if "terms" in stored else ""),
                unsafe_allow_html=True)
