import logging
import re
from datetime import datetime
from typing import Any, Dict, List, Optional, Set

import pandas as pd

from background_gene_set import BackgroundGeneSet
from enrichment import Enrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary

logger = logging.getLogger(__name__)


class IterativeEnrichment:
    """
    Wrapper for iterative gene set enrichment.

    :param gene_set: List of gene identifiers to analyze.
    :type gene_set: GeneSet
    :param background_gene_set: Background genes file (one gene per line).
    :type background_gene_set: BackgroundGeneSet
    :param gene_set_library: Gene set library GMT file.
    :type gene_set_library: GeneSetLibrary
    :param p_value_method_name: Name of p-value calculation method to pass to Enrichment.
    :type p_value_method_name: str
    :param p_threshold: P-value cutoff for including terms.
    :type p_threshold: float
    :param max_iterations: Maximum number of enrichment iterations to run. None means no limit.
    :type max_iterations: Optional[int]
    """

    def __init__(
        self,
        gene_set: GeneSet,
        gene_set_library: GeneSetLibrary,
        background_gene_set: BackgroundGeneSet,
        min_term_size: int = 10,
        max_term_size: int = 600,
        p_value_method_name: str = "Fisher's Exact Test",
        name: str = None,
        p_threshold: float = 0.01,
        max_iterations: Optional[int] = None,
        min_overlap: int = 1,
    ) -> None:
        """
        Initialize iterative enrichment.

        :param gene_set: Input gene set
        :param gene_set_library: Gene set library
        :param background_gene_set: Background gene set
        :param min_term_size: Minimum term size
        :param max_term_size: Maximum term size
        :param p_value_method_name: P-value calculation method
        :param name: Name for the enrichment
        :param p_threshold: P-value cutoff for including terms
        :param max_iterations: Maximum number of iterations (None for no limit)
        :param min_overlap: Minimum overlap size required for terms
        """
        self.gene_set = gene_set
        self.gene_set_library = gene_set_library
        self.min_term_size = min_term_size
        self.max_term_size = max_term_size
        self.background_gene_set = background_gene_set
        self.p_value_method_name: str = p_value_method_name
        self.p_threshold: float = p_threshold
        self.max_iterations: Optional[int] = max_iterations
        self.min_overlap: int = min_overlap
        self.name = (
            name
            if name
            else f"{gene_set.name}_{gene_set_library.name}_{background_gene_set.name}_{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        )
        self._results: List[Dict[str, Any]] = self._compute_enrichment()

    @property
    def results(self) -> List[Dict[str, Any]]:
        """
        The list of iteration records for this enricher.

        :returns: List of dictionaries with keys [iteration, term, library, p-value, genes]
        :rtype: List[Dict[str, Any]]
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
        Perform iterative enrichment, peeling off top terms until no further terms pass p-value threshold.

        :returns: List of iteration records
        :rtype: List[Dict[str, Any]]
        """
        remaining: Set[str] = set(self.gene_set.genes)
        iteration: int = 1
        records: List[Dict[str, Any]] = []

        while True:
            if not remaining:
                logger.info("No genes left; stopping iterative enrichment.")
                break
            if self.max_iterations is not None and iteration > self.max_iterations:
                logger.warning("Reached max_iterations; stopping iterative enrichment.")
                break

            # Create a new GeneSet with the remaining genes, but avoid re-validation
            # since these genes were already validated when the original gene set was created
            current_set = GeneSet(
                list(remaining), 
                set(self.background_gene_set.genes),
                hgcn=False,  # Don't re-validate against background
                format=False  # Don't re-format genes
            )
            try:
                enr = Enrichment(
                    gene_set=current_set,
                    gene_set_library=self.gene_set_library,
                    min_term_size=self.min_term_size,
                    max_term_size=self.max_term_size,
                    background_gene_set=self.background_gene_set,
                    p_value_method_name=self.p_value_method_name,
                )
            except Exception as e:
                logger.error(f"Enrichment failed at iteration {iteration}: {e}")
                break

            results = enr.results
            if not results:
                logger.info("No enrichment results; terminating.")
                break

            # Filter results by minimum overlap (same as regular mode)
            filtered_results = [
                result for result in results
                if (result.get("overlap_size", "").split("/")[0].isdigit() and 
                    int(result.get("overlap_size", "").split("/")[0]) >= self.min_overlap)
            ]
            
            if not filtered_results:
                logger.info(f"No results meet minimum overlap requirement ({self.min_overlap}); terminating.")
                break

            top = filtered_results[0]
            pval = top.get("p-value")
            if pval is None or pval >= self.p_threshold:
                logger.info("Top term p-value >= threshold; terminating.")
                break

            # Get overlap size from overlap_size field (format: "3/50")
            overlap_size_str = top.get("overlap_size", "0/0")
            try:
                overlap_count = int(overlap_size_str.split("/")[0])
            except (ValueError, IndexError):
                logger.warning(f"Could not parse overlap_size: {overlap_size_str}")
                overlap_count = 0
            
            # Get overlap genes for the record
            overlap_data = top.get("overlap", [])
            if isinstance(overlap_data, str):
                # If overlap is a string, try to parse it
                genes_in_term = set(overlap_data.split(',') if overlap_data else [])
            elif isinstance(overlap_data, list):
                genes_in_term = set(overlap_data)
            else:
                logger.warning(f"Unexpected overlap data type: {type(overlap_data)}")
                genes_in_term = set()
            
            # Debug logging
            logger.info(f"Iteration {iteration}: Top term '{top.get('term', '')}' has p-value {pval} and overlap size {overlap_count}")
            logger.info(f"Overlap genes: {genes_in_term}")
            logger.info(f"Minimum overlap requirement: {self.min_overlap}")
            
            # Note: We already filtered by minimum overlap above, so this check is redundant but kept for safety
            if overlap_count < self.min_overlap:
                logger.info(f"Top term overlap size ({overlap_count}) < minimum overlap requirement ({self.min_overlap}); terminating.")
                break
            
            record: Dict[str, Any] = {
                "iteration": iteration,
                "term": top.get("term", ""),
                "library": self.gene_set_library.name,
                "p-value": pval,
                "overlap_size": top.get("overlap_size", "0/0"),
                "genes": sorted(genes_in_term),
            }
            records.append(record)
            remaining -= genes_in_term
            iteration += 1

        return records

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert iteration records to a pandas DataFrame.

        :returns: DataFrame of iteration records
        :rtype: pandas.DataFrame
        """
        return pd.DataFrame(self.results)

    def to_tsv(self) -> str:
        """
        Serialize iteration records as a TSV string.

        :returns: TSV-formatted string
        :rtype: str
        """
        return self.to_dataframe().to_csv(sep="\t", index=False)

    def to_json(self) -> str:
        """
        Serialize iteration records as a JSON-formatted string.

        :returns: JSON-formatted string
        :rtype: str
        """
        import json

        return json.dumps(self.results, indent=2)

    def to_dot(self) -> str:
        """
        Generate a valid Graphviz DOT for the iterative enrichment network,
        with sanitized, quoted IDs, semicolons, and no duplicates.
        """

        def _sanitize_id(raw: str) -> str:
            """
            Convert raw label into a valid DOT node ID: replace non-alphanumeric with underscores.
            Collapse multiple underscores and strip leading/trailing underscores.
            """
            # replace non-word characters with underscore
            s = re.sub(r"\W+", "_", raw)
            # collapse underscores
            s = re.sub(r"_+", "_", s)
            return s.strip("_")

        def _format_term_name(term_name: str) -> str:
            """
            Convert underscores to spaces in term names for better readability.
            """
            return term_name.replace("_", " ")

        nodes: Set[str] = set()
        edges: Set[str] = set()

        # Build nodes and edges
        for rec in self.results:
            term_label = rec.get("term", "")
            # Format term name for display (convert underscores to spaces)
            formatted_term_label = _format_term_name(term_label)
            # sanitize and quote term ID
            raw_id = f"term_{rec['iteration']}_{term_label}"
            term_id = _sanitize_id(raw_id)
            term_node = (
                f'"{term_id}" '
                f'[label="{formatted_term_label}", '
                f'style=filled, fontcolor="white", type="term"];'
            )
            nodes.add(term_node)

            for gene in rec.get("genes", []):
                gene_id = _sanitize_id(f"gene_{gene}")
                gene_node = f'"{gene_id}" [label="{gene}", type="gene"];'
                nodes.add(gene_node)
                # create edge with quoted IDs and semicolon
                edge = f'"{gene_id}" -- "{term_id}";'
                edges.add(edge)

        # Assemble DOT
        lines: List[str] = []
        lines.append("graph iterative_enrichment {")
        lines.append("  graph [layout=neato];")
        lines.append("  node [shape=ellipse];")

        # add nodes
        for node in sorted(nodes):
            lines.append(f"  {node}")
        # add edges
        for edge in sorted(edges):
            lines.append(f"  {edge}")

        lines.append("}")
        return "\n".join(lines)
