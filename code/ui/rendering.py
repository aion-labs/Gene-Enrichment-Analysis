import logging
from math import log10

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from streamlit import session_state as state

from enrichment import Enrichment
from iter_enrichment import IterativeEnrichment
from ui.dot_utils import dot_to_plotly
from ui.utils import download_link

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def render_table(result: pd.DataFrame) -> None:
    """
    Render a styled DataFrame within the Streamlit app.

    This function takes a DataFrame containing results data and applies custom formatting before
    displaying it within the Streamlit app using the `st.dataframe` function. Custom formatting includes
    specific number formatting for 'p-value' and 'fdr' columns.

    :param result: The DataFrame containing results data to display.
    """

    logger.info("Rendering DataFrame in Streamlit app.")

    def custom_format(n):
        if n > 0.001:
            return f"{n:.3f}"
        else:
            return f"{n:.3e}"

    # Create a copy of the dataframe to avoid modifying the original
    display_df = result.copy()
    
    # Convert underscores to spaces in term names for better readability
    if 'term' in display_df.columns and not display_df.empty:
        # Ensure term column contains strings and handle NaN/None values
        display_df['term'] = display_df['term'].fillna('').astype(str).str.replace('_', ' ', regex=False)

    st.dataframe(
        display_df.style.format({"p-value": custom_format, "fdr": custom_format}),
        use_container_width=True,
        column_config={
            "rank": None,
            "term": "Term",
            "overlap_size": "Overlap size",
            "overlap": "Overlap (click to expand)",
            "description": None,
            "p-value": "P-value",
            "fdr": "FDR",
        },
    )


def render_barchart(result: pd.DataFrame, file_name: str = "") -> None:
    """
    Render a bar chart visualization of the results data within the Streamlit app.

    This function takes a DataFrame containing results data, specifically terms and p-values, and creates
    a bar chart using Plotly Express, which is then displayed in the Streamlit app using the `st.plotly_chart` function.

    :param result: The DataFrame containing results data to visualize.
    :param file_name: Optional file name to create unique chart keys.
    """
    logger.info("Rendering bar chart in Streamlit app.")
    bar = result[["term", "p-value"]].copy()
    
    # Convert underscores to spaces in term names for better readability
    if not bar.empty:
        bar['term'] = bar['term'].fillna('').astype(str).str.replace('_', ' ', regex=False)
    
    bar.loc[:, "p-value"] = bar.loc[:, "p-value"].apply(lambda x: -1 * log10(x))
    bar = bar.sort_values(by=["p-value"])
    bar.columns = ["term", "-log10(p-value)"]
    fig = px.bar(
        bar,
        x="-log10(p-value)",
        y="term",
        orientation="h",
        labels={"-log10(p-value)": "−log₁₀(p‐value)", "term": "Term"},
    )
    st.plotly_chart(fig, use_container_width=True, key=f"barchart_{file_name}")


def render_results(result: Enrichment, file_name: str, n_results: int = 10) -> None:
    """
    Render a results section within the Streamlit app.

    This function processes and visualizes the result data within the Streamlit app. It provides
    a table visualization and a bar chart visualization of the top results. Additionally, it
    allows the user to download the results in various formats.

    :param result: The DataFrame containing results data to display.
    :param file_name: The name of the file to be used for downloading results.
    :param n_results: Numbers of results to display
    """
    logger.info(f"Rendering results for file: {file_name}")
    result_df = result.to_dataframe().head(n_results)
    result_df = result_df.set_index("rank")
    st.divider()
    st.subheader(file_name)
    table, bar = st.tabs(["Results", "Bar chart"])
    with table:
        render_table(result_df)
    with bar:
        render_barchart(result_df, file_name)

    st.markdown(
        f'Download results as {download_link(result.to_tsv(), file_name, "tsv")}, {download_link(result.to_json(), file_name, "json")}',
        unsafe_allow_html=True,
    )


def render_validation() -> None:
    """
    Validate and render the input gene set information.

    This function checks the `gene_set` in the session state for duplicates and invalid entries.
    It then provides a feedback to the user in the Streamlit app on the validation results.
    """
    logger.info("Validating and rendering the input gene set information.")
    if "gene_set" in state:
        total = state.gene_set.size
        dups = len(state.gene_set.validation["duplicates"])
        non_gene = len(state.gene_set.validation["non_valid"])
        if dups:
            dups_st = f", ⚠️ {dups} duplicates"
        else:
            dups_st = ""
        if non_gene:
            non_gene_st = f", ⛔ {non_gene} non valid"
        else:
            non_gene_st = ""
        caption = f"{total} genes{dups_st}{non_gene_st}"
        with st.expander(caption):
            if dups:
                st.data_editor(
                    pd.json_normalize(state.gene_set.validation)["duplicates"],
                    column_config={
                        "duplicates": st.column_config.ListColumn(
                            dups_st[2:], width="large"
                        ),
                    },
                    hide_index=True,
                )
            if non_gene:
                st.data_editor(
                    pd.json_normalize(state.gene_set.validation)["non_valid"],
                    column_config={
                        "non_valid": st.column_config.ListColumn(
                            non_gene_st[2:], width="large"
                        ),
                    },
                    hide_index=True,
                )


def render_iter_table(result: pd.DataFrame) -> None:
    """
    Render a styled DataFrame of iterative enrichment results.

    :param result: DataFrame indexed by 'iteration' with columns ['term', 'p-value', 'overlap_size', 'genes'].
    :type result: pandas.DataFrame
    """
    logger.info("Rendering iterative results table.")

    # Create a copy to avoid modifying the original
    df = result.copy()
    
    # Convert underscores to spaces in term names for better readability
    if 'term' in df.columns and not df.empty:
        df['term'] = df['term'].fillna('').astype(str).str.replace('_', ' ', regex=False)

    # Rename 'genes' column for clarity
    df = df.rename(columns={"genes": "Genes removed"})

    # Apply custom formatting to p-value column
    def custom_format(n):
        if n > 0.001:
            return f"{n:.3f}"
        return f"{n:.2e}"

    styled = df.style.format({"p-value": custom_format})

    st.dataframe(
        styled,
        use_container_width=True,
        column_config={
            "term": "Term",
            "p-value": "P-value",
            "overlap_size": "Overlap size",
            "Genes removed": "Genes removed",
        },
    )


def render_iter_barchart(result: pd.DataFrame, file_name: str = "") -> None:
    """
    Render a bar chart visualization of iterative enrichment results.

    :param result: DataFrame indexed by 'iteration' with at least a 'p-value' column.
    :type result: pandas.DataFrame
    :param file_name: Optional file name to create unique chart keys.
    :type file_name: str
    """
    logger.info("Rendering iterative bar chart.")

    # Prepare bar plot data
    bar = result.reset_index()[["iteration", "term", "p-value"]].copy()
    
    # Convert underscores to spaces in term names for better readability
    if not bar.empty:
        bar['term'] = bar['term'].fillna('').astype(str).str.replace('_', ' ', regex=False)
    
    bar["-log10(p-value)"] = bar["p-value"].apply(
        lambda x: -log10(x) if x and x > 0 else None
    )
    bar = bar.sort_values(by=["p-value"], ascending=False)

    fig = px.bar(
        bar,
        x="-log10(p-value)",
        y="term",
        orientation="h",
        hover_data=["p-value"],
        labels={"term": "Term", "-log10(p-value)": "-log10(p-value)"},
        title="Iterative Enrichment p-value per Iteration",
    )
    st.plotly_chart(fig, use_container_width=True, key=f"iter_barchart_{file_name}")


def render_iter_results(result: IterativeEnrichment, file_name: str) -> None:
    """
    Render a results section for iterative enrichment within the Streamlit app.

    :param result: IterativeEnrichment object containing iteration records.
    :type result: src.iter_enrichment.IterativeEnrichment
    :param file_name: Name of the library or section header.
    :type file_name: str
    """
    logger.info(f"Rendering iterative results for {file_name}.")

    st.divider()
    st.subheader(file_name)

    # Tabs for table and chart
    table_tab, bar_tab = st.tabs(["Iterations", "Bar chart"])
    df = result.to_dataframe()

    if df.empty or "iteration" not in df.columns:
        logger.warning("No iterative enrichment results.")
        st.warning("No iterative enrichment results.")
    else:
        df = df.set_index("iteration")
        with table_tab:
            render_iter_table(df)

        with bar_tab:
            render_iter_barchart(df, file_name)

        # Download links
        st.markdown(
            f'Download iterative results as {download_link(result.to_tsv(), file_name, "tsv")}, '
            f'{download_link(result.to_json(), file_name, "json")}',
            unsafe_allow_html=True,
        )


def render_network(dot: str, title: str = "Iterative Enrichment Network") -> None:
    """
    Render a bipartite network graph of iterative enrichment across libraries.

    :param dot: Graphviz DOT-format string.
    :type dot: str
    :param title: Header title for the network graph.
    :type title: str
    """
    logger.info("Rendering iterative network graph.")

    st.divider()
    st.subheader(title)
    # st.graphviz_chart(dot, use_container_width=True)

    st.plotly_chart(dot_to_plotly(dot), use_container_width=True, key="network_chart")
    # Offer DOT download and AI analysis prompt download
    st.markdown(
        f'Download network graph as {download_link(dot, "iterative_network", "dot")}, '
        f'Download {download_link(generate_ai_analysis_prompt(dot), "ai_analysis_prompt", "txt")} for AI analysis',
        unsafe_allow_html=True,
    )


def generate_ai_analysis_prompt(dot_content: str) -> str:
    """
    Generate an AI analysis prompt based on the DOT network content.
    
    :param dot_content: The DOT network content
    :return: Formatted AI analysis prompt
    """
    prompt = """You are a computational biologist analyzing an iterative gene set enrichment network. Iterative means that for each library source, the top enriched gene set is found, saved and the genes are removed from the initial gene list. The remaining genes are then tested for enrichment again. Thus, each gene appears only once per library source tested but can be linked to multiple terms that originate from different libraries. The names of the libraries are indicated in the legend. The results are provided as a DOT network represents the relationships between genes and gene sets.

**NETWORK STRUCTURE:**
"""
    
    # Add the DOT content
    prompt += dot_content + "\n\n"
    
    # Add the rest of the prompt
    prompt += """**NETWORK INTERPRETATION GUIDE:**

**Node Types:**
- **Gene nodes** (type="gene"): Individual genes from the input gene set
- **Term nodes** (type="term"): Enriched biological pathways/processes/gene sets
- **Colors**: Term nodes are colored by library source

**Edge Connections:**
- Each edge (gene -- term) represents a gene's membership in that biological term
- Thicker/more connections indicate genes that are part of multiple enriched processes
- The network topology shows which genes are central to the biological response

**Iteration Information:**
- Term nodes are labeled with iteration numbers (term_1_, term_2_, etc.)
- Earlier iterations (lower numbers) represent the most significant enrichments
- Ignore iteration number for biological interpretation, instead focus on network hubs densly connected for interesting biological insights

**ANALYSIS REQUEST:**

Please analyze this network and provide:

1. **Key Biological Insights:**
   - What are the most significant biological processes/pathways identified?
   - Which genes appear to be central to the biological response?
   - What patterns emerge from the iteration sequence?

2. **Network Topology Analysis:**
   - Which genes are "hub" genes (connected to multiple terms)?
   - Are there distinct clusters or modules in the network?
   - What does the connectivity pattern suggest about biological organization?

3. **Biological Hypothesis:**
   - Based on the network structure, what biological hypothesis can you generate?
   
4. **Estimated Experimental Context**
   - Can you hypothesize on what was the experiment that generated this result?

**RESPONSE STRUCTURE:**
- **Executive Summary** (2-3 sentences)
- **Key Biological Processes Identified**
- **Network Topology Insights**
- **Proposed Biological Hypothesis**
- **Estimated Experimental Context**

**IMPORTANT NOTES:**
- Focus on biological interpretation, not statistical significance or iteration number
- Consider the functional relationships between connected genes/terms
- Look for unexpected connections that might reveal novel biology
- Consider the broader biological context and literature"""
    
    return prompt
