import json
import logging
import multiprocessing as mp
from datetime import datetime
from typing import Any, Dict, List, Tuple

import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def compute_pvalue(
    args: Tuple[GeneSet, BackgroundGeneSet, dict, str, int],
) -> Tuple[str, str, str, List[str], float]:
    """
    Computes the p-value for a given term using Fisher's exact test.
    This function is intended to be used with multiprocessing.Pool.map(),
    which requires functions to take a single argument. Therefore, the
    inputs are passed as a single tuple.

    Args:
        args: A tuple containing the following elements:
            - gene_set (GeneSet): The input gene set
            - background_gene_set (BackgroundGeneSet): The background gene set
            - term (dict): A dictionary representing a term in the gene set library
            - p-value method (str): The name of the method to calculate p-value
            - library_background_size (int): The size of the intersection between background and library genes

    Returns:
        A tuple containing the following elements:
            - term name (str): The name of the term
            - overlap size (int): The size of the overlap between the gene set and the term
            - term description (str): The description of the term
            - overlap genes (list): A list of genes in the overlap
            - p_value (float): The p-value computed by Fisher's exact test
    """
    gene_set, background_gene_set, term, p_value_method_name, library_background_size = args
    term_genes = set(term["genes"])
    n_term_genes = len(term_genes)
    overlap = gene_set.genes & term_genes
    n_overlap = len(overlap)

    # Build contingency table for Fisher's exact test using library-specific background
    contingency_table = [
        [n_overlap, n_term_genes - n_overlap],
        [
            gene_set.size - n_overlap,
            library_background_size - n_term_genes - gene_set.size + n_overlap,
        ],
    ]

    if p_value_method_name == "Fisher's Exact Test":
        _, p_value = fisher_exact(contingency_table)
    elif p_value_method_name == "Chi-squared Test":
        chi2, p_value, _, _ = chi2_contingency(contingency_table)
    elif p_value_method_name == "Hypergeometric Test":
        p_value = hypergeom.sf(
            n_overlap - 1, library_background_size, n_term_genes, gene_set.size
        )
    else:
        logger.error(f"Unsupported p_value_method: {p_value_method_name}")
        raise ValueError(f"Unsupported p_value_method: {p_value_method_name}")
    return (
        term["name"],
        f'{len(overlap)}/{len(term["genes"])}',
        term["description"],
        sorted(list(overlap)),
        p_value,
    )


class Enrichment:
    """
    Class for gene set enrichment analysis results.
    """

    def __init__(
        self,
        gene_set: GeneSet,
        gene_set_library: GeneSetLibrary,
        background_gene_set: BackgroundGeneSet,
        min_term_size: int = 10,
        max_term_size: int = 1000,
        p_value_method_name="Fisher's Exact Test",
        name: str = None,
    ):
        """
        Initialize the class with gene set, gene set library, and background gene set.

        Args:
            gene_set: Input gene set
            gene_set_library: Gene set library
            background_gene_set: Background gene set
        """
        self.gene_set = gene_set
        self.gene_set_library = gene_set_library
        self.min_term_size = min_term_size
        self.max_term_size = max_term_size
        self.background_gene_set = background_gene_set
        self.p_value_method_name = p_value_method_name
        self.name = (
            name
            if name
            else f"{gene_set.name}_{gene_set_library.name}_{min_term_size}-{max_term_size}_{background_gene_set.name}_{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        )
        self._results: List[Dict[str, Any]] = self._compute_enrichment()

    @property
    def results(self) -> List[Dict[str, Any]]:
        """
        Getter for _results.

        Returns:
            A list containing dictionaries of enrichment results
        """
        return self._results

    @results.setter
    def results(self, value: List[Dict[str, Any]]) -> None:
        """
        Setter for _results.

        Args:
            value: A list containing dictionaries of enrichment results
        """
        self._results = value

    def _compute_enrichment(self) -> List[Dict[str, Any]]:
        """
        Computes gene set enrichment analysis.

        Returns:
            A list containing dictionaries of enrichment results
        """
        results = []
        logger.info(f"Calculating p-values for {self.gene_set_library.name}")
        cpu_count = mp.cpu_count() - 2
        parallel_results = []  # Initialize outside try block
        
        with mp.Pool(cpu_count) as pool:
            logger.info(f"Initializing the MP pool with {cpu_count} CPUs")
            try:
                # Calculate the size of the intersection between the background gene set and the library's unique genes
                library_background_size = len(self.background_gene_set.genes & self.gene_set_library.unique_genes)
                logger.info(f"Library-specific background size: {library_background_size} genes (intersection of {self.background_gene_set.size} background genes and {self.gene_set_library.size} library genes)")
                parallel_results = pool.map(
                    compute_pvalue,
                    [
                        (
                            self.gene_set,
                            self.background_gene_set,
                            term,
                            self.p_value_method_name,
                            library_background_size,
                        )
                        for term in self.gene_set_library.library if self.min_term_size <= term["size"] <= self.max_term_size
                    ],
                )
            except Exception as e:
                logging.exception("An error occurred: %s", e)
                return results  # Return empty results if computation failed
            finally:
                pool.close()
                pool.join()
                logger.info(f"Releasing {cpu_count} CPUs from the MP pool")

        # Check if we have results to process
        if not parallel_results:
            logger.warning(f"No results obtained for {self.gene_set_library.name}")
            logger.info(f"Library has {len(self.gene_set_library.library)} total terms")
            logger.info(f"Terms within size range [{self.min_term_size}, {self.max_term_size}]: {len([t for t in self.gene_set_library.library if self.min_term_size <= t['size'] <= self.max_term_size])}")
            logger.info(f"Input gene set size: {self.gene_set.size}")
            logger.info(f"Background gene set size: {self.background_gene_set.size}")
            return results

        # Separate results and p_values for convenience
        p_values = [result[-1] for result in parallel_results]
        # Adjust p-values for multiple testing
        _, p_values_adjusted, _, _ = multipletests(p_values, method="fdr_bh")
        # Rank terms based on their p-values
        ranked_terms = sorted(
            list(enumerate(parallel_results)), key=lambda x: p_values[x[0]]
        )

        # Format results into a sorted list
        for i, result in ranked_terms:
            term_name, overlap_size, term_description, overlap_genes, _ = result
            results.append(
                {
                    "term": term_name,
                    "rank": i + 1,
                    "description": term_description,
                    "overlap": overlap_genes,
                    "overlap_size": overlap_size,
                    "p-value": p_values[i],
                    "fdr": p_values_adjusted[i],
                }
            )
        return results

    def to_dataframe(self):
        """Return the enrichment results as a pandas dataframe."""
        return pd.DataFrame(
            {
                "rank": [result["rank"] for result in self.results],
                "term": [result["term"] for result in self.results],
                "overlap": [result["overlap"] for result in self.results],
                "overlap_size": [result["overlap_size"] for result in self.results],
                "p-value": [result["p-value"] for result in self.results],
                "fdr": [result["fdr"] for result in self.results],
            }
        )

    def to_json(self):
        """Return the enrichment results as a JSON string."""
        return json.dumps(self.results, indent=4, separators=(",", ": "))

    def to_html(self):
        """Return the enrichment results as an HTML page."""
        return self.to_dataframe().to_html()

    def to_tsv(self):
        """Return the enrichment results as a TSV spreadsheet."""
        return self.to_dataframe().to_csv(sep="\t")

    def to_snapshot(self) -> Dict:
        """Return the snapshot of input parameters and the enrichment results as a JSON string."""
        # Calculate library-specific background size for the snapshot
        library_background_size = len(self.background_gene_set.genes & self.gene_set_library.unique_genes)
        return {
            "input_gene_set": list(self.gene_set.genes),
            "background": self.background_gene_set.name,
            "background_size": self.background_gene_set.size,
            "library_background_size": library_background_size,
            "library_size": self.gene_set_library.size,
            self.gene_set_library.name: self.results,
        }
