import json
import logging
import math
from io import StringIO
from pathlib import Path
from typing import Dict, List, Set

import streamlit as st
from PIL import Image
from streamlit import session_state as state

# Existing imports
from background_gene_set import BackgroundGeneSet
from enrichment import Enrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from iter_enrichment import IterativeEnrichment
from ui.dot_utils import merge_iterative_dot
from ui.helpers import input_example, update_text_widgets, convert_and_validate_gene_input, display_conversion_results
from ui.processing import collect_results
from ui.rendering import (
    render_iter_results,
    render_network,
    render_results,
    render_validation,
)
from ui.utils import download_link, update_aliases

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parent.parent

st.set_page_config(
    page_title="Enrichment Analysis", layout="wide", initial_sidebar_state="expanded"
)


def _ensure_base_state():
    if "enrich" not in state:
        state.enrich = {}
    if "iter_results" not in state:
        state.iter_results: Dict[str, List[dict]] = {}
    if "iter_graph_parts" not in state:
        state.iter_graph_parts: Dict[str, dict] = {}
    if "results_ready" not in state:
        state.results_ready = False
    if "iter_min_overlap" not in state:
        state.iter_min_overlap = 3
    if "min_term_size" not in state:
        state.min_term_size = 10
    if "max_term_size" not in state:
        state.max_term_size = 1000
    if "iter_ready" not in state:
        state.iter_ready = False
    if "selected_dot_paths" not in state:
        state.selected_dot_paths = []
    if "network_generated" not in state:
        state.network_generated = False  # flag to prevent clearing after checkbox changes


def _build_iterative_tables_download(all_iter_results: Dict[str, List[dict]]) -> str:
    rows = ["library\titeration\tterm\t p-value\t-log10_p\tgenes"]
    for lib, records in all_iter_results.items():
        for rec in records:
            p = rec.get("p-value", float("nan"))
            rows.append(
                f"{lib}\t{rec['iteration']}\t{rec['term']}\t{p}\t"
                f"{(-math.log10(p) if p and p>0 else '')}\t{','.join(rec.get('genes', []))}"
            )
    return "\n".join(rows)


def main() -> None:
    logger.info("Starting the Streamlit app")
    st.sidebar.image(
        Image.open(ROOT / "code" / "static" / "logo.png"),
        caption="Iterative Enrichment Analysis",
    )
    st.sidebar.title("Enrichment analysis")
    st.sidebar.write(
        """This app tests the input gene list against selected pathway libraries and iteratively removes genes from the top hit at each step. The process stops when no term passes the p-value cutoff. Results include a ranked table, bar chart, and network graph that highlight both primary and secondary signals."""
    )

    _ensure_base_state()

    mode = st.radio(
        "Mode",
        ["Regular", "Iterative"],
        index=1,
        horizontal=True,
        key="analysis_mode",
    )
    st.subheader(f"Enrichment analysis — {mode} mode")

    state.lib_mapper = update_aliases("libraries")
    state.bg_mapper = update_aliases("backgrounds")
    state.advanced_settings_changed = False
    state.bt_submit_disabled = True

    analysis, advanced = st.tabs(["Analysis", "Advanced settings"])

    with analysis:
        col_input, col_settings = st.columns([5, 7])
        with col_input:
            # Gene input format selector
            if 'gene_input_format' not in state:
                state.gene_input_format = 'symbols'
            
            format_col1, format_col2 = st.columns([1, 3])
            with format_col1:
                state.gene_input_format = st.selectbox(
                    "Input Format",
                    ["symbols", "entrez_ids"],
                    index=0 if state.gene_input_format == 'symbols' else 1,
                    format_func=lambda x: "Gene Symbols" if x == "symbols" else "Entrez IDs"
                )
            
            with format_col2:
                if state.gene_input_format == 'symbols':
                    st.caption("💡 **Gene Symbols:** Enter official gene symbols (e.g., TP53, BRCA1)")
                else:
                    st.caption("💡 **Entrez IDs:** Enter numeric Entrez IDs (e.g., 7157, 672)")
            
            # Gene input area
            placeholder_text = (
                "Enter gene symbols (e.g., TP53, BRCA1) - one per line" 
                if state.gene_input_format == 'symbols' 
                else "Enter Entrez IDs (e.g., 7157, 672) - one per line"
            )
            
            st.text_area(
                "Input a set of genes",
                key="gene_set_input",
                height=400,
                placeholder=placeholder_text,
                label_visibility="collapsed",
            )
            st.text_input(
                "Gene set name",
                key="gene_set_name",
                placeholder="Input a gene set name",
                label_visibility="collapsed",
            )
            gene_files = [
                str(f).replace(f"{ROOT}/data/gene_lists/", "")
                for f in (ROOT / "data" / "gene_lists").rglob("*.txt")
            ]
            st.selectbox(
                "Or select a file from the `data` folder",
                ["Select ..."] + gene_files,
                index=0,
                on_change=update_text_widgets,
                key="selected_file",
            )
        with col_settings:
            state.background_set = st.selectbox(
                "Background gene list", state.bg_mapper.keys()
            )
            st.caption("Specifies the background list of genes...")
            state.libraries = st.multiselect(
                "Select libraries",
                state.lib_mapper.keys(),
                # default=list(state.lib_mapper.keys())
            )
            if state.libraries:
                state.gene_set_libraries = [
                    GeneSetLibrary(
                        str(ROOT / "data" / "libraries" / state.lib_mapper[lib]),
                        name=lib,
                    )
                    for lib in state.libraries
                ]
                # filter out oversized terms by max_term_size setting
                for gsl in state.gene_set_libraries:
                    filtered_terms = [
                        t for t in gsl.library if t["size"] <= state.iter_max_term_size
                    ]
                    gsl.library = filtered_terms
                    gsl.num_terms = len(filtered_terms)
                    gsl.unique_genes = gsl.compute_unique_genes()
                    gsl.size = len(gsl.unique_genes)
            else:
                state.gene_set_libraries = []
            if state.background_set:
                state.background_gene_set = BackgroundGeneSet(
                    str(
                        ROOT
                        / "data"
                        / "backgrounds"
                        / state.bg_mapper[state.background_set]
                    )
                )
                if state.gene_set_input:
                    # Convert and validate gene input based on selected format
                    converted_symbols, unrecognized_entrez, unrecognized_symbols, stats = convert_and_validate_gene_input(
                        state.gene_set_input, 
                        state.gene_input_format
                    )
                    
                    # Display conversion results
                    display_conversion_results(converted_symbols, unrecognized_entrez, unrecognized_symbols, stats, state.gene_input_format)
                    
                    # Create gene set with converted symbols
                    if converted_symbols:
                        state.gene_set = GeneSet(
                            converted_symbols,
                            state.background_gene_set.genes,
                            state.gene_set_name,
                        )
                    else:
                        state.gene_set = None
            if mode == "Iterative":
                st.markdown("**Iterative parameters**")
                state.iter_p_threshold = st.number_input(
                    "P-value threshold",
                    min_value=1e-10,
                    max_value=0.5,
                    value=0.01,
                    step=0.001,
                    format="%.4f",
                )
                state.iter_max_iter = st.number_input(
                    "Max iterations (0 = no limit)",
                    min_value=0,
                    max_value=500,
                    value=10,
                    step=1,
                )
                state.iter_min_overlap = st.number_input(
                    "Minimum overlap with gene set",
                    min_value=1,
                    value=3,
                    step=1,
                )
                state.iter_min_term_size, state.iter_max_term_size = st.slider(
                    "Minimum and maximum term size",
                    min_value=1,
                    value=(10, 1000),
                    step=10,
                    max_value=5000
                )


    col_sub, col_example, _ = st.columns([9, 8, 29])
    ready_common = all(
        getattr(state, k, None)
        for k in ["gene_set", "background_gene_set", "gene_set_libraries"]
    )
    with col_sub:
        if mode == "Regular":
            state.bt_submit_disabled = not ready_common
            bt_submit = st.button(
                "Submit", disabled=state.bt_submit_disabled, key="bt_reg"
            )
        else:
            state.bt_iter_disabled = not ready_common
            bt_iter = st.button(
                "Submit",
                disabled=state.bt_iter_disabled,
                key="bt_iter",
            )
    with col_example:
        st.button("Input an example", on_click=input_example)

    with advanced:
        if mode == "Regular":
            n_results = st.slider(
                "Number of results to display", 1, 100, 10, 1, key="n_res"
            )
        else:
            st.slider(
                "Number of results to display (regular only)",
                1,
                100,
                10,
                1,
                disabled=True,
            )
        # Use widget key to set session_state; do not assign to state directly
        st.selectbox(
            "P-value calculation method",
            ["Fisher's Exact Test", "Hypergeometric Test", "Chi-squared Test"],
            key="p_val_method",
        )
        if state.p_val_method != "Fisher's Exact Test":
            state.advanced_settings_changed = True
        # Background gene list format selector
        bg_format_col1, bg_format_col2 = st.columns([1, 3])
        with bg_format_col1:
            if 'bg_input_format' not in state:
                state.bg_input_format = 'symbols'
            
            state.bg_input_format = st.selectbox(
                "Background Format",
                ["symbols", "entrez_ids"],
                index=0 if state.bg_input_format == 'symbols' else 1,
                format_func=lambda x: "Gene Symbols" if x == "symbols" else "Entrez IDs"
            )
        
        with bg_format_col2:
            if state.bg_input_format == 'symbols':
                st.caption("💡 **Background Symbols:** Upload file with gene symbols")
            else:
                st.caption("💡 **Background Entrez IDs:** Upload file with Entrez IDs")
        
        state.bg_custom = st.file_uploader(
            "Upload your background gene list", type=[".txt"]
        )
        if state.bg_custom:
            bgf = (ROOT / "data" / "backgrounds" / state.bg_custom.name).open("wb")
            bgf.write(state.bg_custom.getvalue())
            state.advanced_settings_changed = True
            
            # Create background gene set with the uploaded file and selected format
            try:
                state.background_gene_set = BackgroundGeneSet(
                    str(ROOT / "data" / "backgrounds" / state.bg_custom.name),
                    name=state.bg_custom.name.replace(".txt", ""),
                    input_format=state.bg_input_format
                )
                st.success(f"✅ Background gene list loaded: {state.background_gene_set.size} genes")
            except Exception as e:
                st.error(f"❌ Error loading background gene list: {str(e)}")
        state.libs_custom = st.file_uploader(
            "Upload gene set libraries",
            type=[".gmt"],
            accept_multiple_files=True,
            on_change=update_aliases,
            args=("libraries",),
        )
        if state.libs_custom:
            for libf in state.libs_custom:
                lf = (ROOT / "data" / "libraries" / libf.name).open("wb")
                lf.write(libf.getvalue())
                state.advanced_settings_changed = True
        if state.advanced_settings_changed:
            if st.button("Apply settings"):
                logger.info("Applied custom settings")
                st.success("Settings applied")
        else:
            with st.empty():
                st.button("Apply settings", disabled=True)

    # Regular execution
    if mode == "Regular" and "bt_submit" in locals() and bt_submit:
        logger.info("Running regular enrichment")
        render_validation()
        if state.gene_set_input and ready_common:
            n_genes = len(state.gene_set_input.split())
            if n_genes <= 100 or n_genes >= 5000:
                warn = "small" if n_genes <= 100 else "big"
                s = "s" if str(n_genes)[-1] != "1" else ""
                st.warning(f"You've entered {n_genes} gene{s}, which may be {warn}...")
            with st.spinner("Calculating enrichment"):
                for gsl in state.gene_set_libraries:
                    enrich = Enrichment(
                        state.gene_set,
                        gsl,
                        state.background_gene_set,
                        state.p_val_method,
                    )
                    state.enrich[gsl.name] = enrich
                    with (ROOT / "results" / f"{enrich.name}.json").open("w") as js:
                        json.dump(enrich.to_snapshot(), js)
                state.results_ready = True
        else:
            if not state.gene_set_input:
                st.error("Please input genes")
            if not state.gene_set_libraries:
                st.error("No libraries selected")
            if not getattr(state, "background_gene_set", None):
                st.error("No background gene list selected")
    if mode == "Regular" and state.results_ready:
        logger.info("Displaying regular results")
        st.markdown(
            f"Download all results as {download_link(collect_results(state.enrich), 'results','tsv')}",
            unsafe_allow_html=True,
        )
        for lib in state.enrich:
            render_results(state.enrich[lib], lib, n_results)
        state.results_ready = False

    # Iterative execution
    if mode == "Iterative" and "bt_iter" in locals() and bt_iter:
        logger.info("Running iterative enrichment")
        render_validation()
        if not state.gene_set_input:
            st.error("Please input genes")
        if not state.gene_set_libraries:
            st.error("No libraries selected")
        if not getattr(state, "background_gene_set", None):
            st.error("No background gene list selected")
        if ready_common and state.gene_set_input:
            # clear prior iterative objects/results
            state.iter_enrich = {}
            state.iter_results.clear()
            state.iter_dot = {}

            # --- Inside your enrichment block in streamlit_app.py ---
            state.iter_enrich = {}
            state.iter_dot.clear()

            with st.spinner("Running iterative enrichment"):
                for gsl in state.gene_set_libraries:
                    it = IterativeEnrichment(
                        gene_set=state.gene_set,
                        gene_set_library=gsl,
                        background_gene_set=state.background_gene_set,
                        min_term_size=state.min_term_size,
                        max_term_size=state.max_term_size,
                        p_value_method_name=state.p_val_method,
                        p_threshold=state.iter_p_threshold,
                        max_iterations=(
                            None if state.iter_max_iter == 0 else state.iter_max_iter
                        ),
                    )
                    # store enrichment object and results
                    state.iter_enrich[gsl.name] = it
                    state.iter_results[gsl.name] = it.results
                    state.iter_dot[gsl.name] = it.to_dot()
                    state.iter_results[gsl.name] = [
                        rec
                        for rec in state.iter_results[gsl.name]
                        if len(rec.get("genes", [])) >= state.iter_min_overlap
                    ]

            state.iter_ready = True

    # Iterative rendering
    if mode == "Iterative" and state.iter_ready:
        combined = _build_iterative_tables_download(state.iter_results)
        st.markdown(
            f"Download iterative results as {download_link(combined, 'iterative_results', 'tsv')}",
            unsafe_allow_html=True,
        )

        # callback to keep checkbox state in session
        def toggle_library(lib_name):
            if state[f"use_{lib_name}_in_network"]:
                if lib_name not in state.selected_dot_paths:
                    state.selected_dot_paths.append(lib_name)
            else:
                if lib_name in state.selected_dot_paths:
                    state.selected_dot_paths.remove(lib_name)

        # render each library's results with a persistent checkbox
        for lib, it in state.iter_enrich.items():
            # Monkey-patch filtered results
            it.results = [
                rec for rec in it.results
                if len(rec.get("genes", [])) >= state.iter_min_overlap
            ]
            render_iter_results(it, lib)
            state.setdefault(f"use_{lib}_in_network", False)
            st.checkbox(
                "Use results in network",
                key=f"use_{lib}_in_network",
                on_change=toggle_library,
                args=(lib,),
            )

        # Network section
        st.markdown("---")
        st.header("Network")
        if state.selected_dot_paths:
            st.write("**Selected libraries:**")
            for sel in state.selected_dot_paths:
                st.write(f"- {sel}")
        else:
            st.write("No libraries selected for network generation.")

        # generate or re-display merged network
        if st.button("Generate Network"):
            state.network_generated = True
            selected_dots = {lib: state.iter_dot[lib] for lib in state.selected_dot_paths}
            state.last_merged_dot = merge_iterative_dot(selected_dots)
            render_network(state.last_merged_dot)
        elif state.network_generated:
            render_network(state.last_merged_dot)

    logger.info("Finishing the Streamlit app")


if __name__ == "__main__":
    main()
