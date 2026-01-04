from pathlib import Path
from typing import List, Optional, Dict
import json
import logging
from datetime import datetime

import typer
import pandas as pd

from background_gene_set import BackgroundGeneSet
from enrichment import Enrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from iter_enrichment import IterativeEnrichment
from ui.utils import get_library_name_from_alias

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
app = typer.Typer(
    help="Gene Set Enrichment Analysis CLI",
    add_completion=False,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]}
)


def run_enrichment(
    gene_sets: List[Path],
    background: Path,
    libraries: List[Path],
    gene_format: str,
    p_value_method: str,
    mode: str,
    p_threshold: float,
    min_overlap: int,
    min_term_size: int,
    max_term_size: int,
    max_iterations: int,
    output_dir: Path,
    benchmark: bool = False,
):
    """Run enrichment analysis for the given parameters."""
    
    # Ensure output directory exists (only if provided)
    if output_dir is not None:
        output_dir.mkdir(exist_ok=True)
    
    # Load background gene set
    logger.info(f"Loading background gene set: {background}")
    background_gene_set = BackgroundGeneSet(str(background))
    
    # Load gene set libraries
    logger.info(f"Loading {len(libraries)} gene set libraries")
    gene_set_libraries = []
    for lib_path in libraries:
        # Get library name from alias.json
        lib_name = get_library_name_from_alias(lib_path)
        library = GeneSetLibrary(str(lib_path), name=lib_name)
        # Filter terms by size (same as Streamlit app)
        filtered_terms = [
            t for t in library.library if min_term_size <= t["size"] <= max_term_size
        ]
        library.library = filtered_terms
        library.num_terms = len(filtered_terms)
        library.unique_genes = library.compute_unique_genes()
        library.size = len(library.unique_genes)
        gene_set_libraries.append(library)
        logger.info(f"Loaded library {lib_name}: {library.num_terms} terms")
    
    # Process each gene set
    for gene_set_path in gene_sets:
        logger.info(f"Processing gene set: {gene_set_path}")
        
        # Load gene set
        gene_set_name = gene_set_path.stem
        with open(gene_set_path, 'r') as f:
            gene_input = [line.strip() for line in f if line.strip()]
        
        # Handle gene format conversion
        if gene_format == "entrez_ids":
            # Convert Entrez IDs to gene symbols
            from gene_converter import GeneConverter
            converter = GeneConverter()
            gene_symbols = []
            unrecognized_entrez = []
            
            for gene_id in gene_input:
                symbol = converter.get_symbol(gene_id)
                if symbol:
                    gene_symbols.append(symbol)
                else:
                    unrecognized_entrez.append(gene_id)
            
            logger.info(f"Converted {len(gene_symbols)} Entrez IDs to gene symbols")
            if unrecognized_entrez:
                logger.warning(f"Warning: {len(unrecognized_entrez)} Entrez IDs not found in database")
        else:
            # Use gene symbols directly
            from gene_converter import GeneConverter
            converter = GeneConverter()
            gene_symbols = []
            unrecognized_symbols = []
            
            for gene_id in gene_input:
                mapped_symbol = converter.validate_and_map_symbol(gene_id)
                if mapped_symbol:
                    gene_symbols.append(mapped_symbol)
                else:
                    unrecognized_symbols.append(gene_id)
            
            if unrecognized_symbols:
                logger.warning(f"Warning: {len(unrecognized_symbols)} gene symbols not found in database")
        
        # Create GeneSet object
        gene_set = GeneSet(
            gene_symbols,
            background_gene_set.genes,
            gene_set_name,
            hgcn=False,
            format=False,
        )
        
        logger.info(f"Gene set {gene_set_name}: {gene_set.size} genes")
        
        # Check gene list size limit (800 genes maximum)
        if gene_set.size > 800:
            logger.error(f"❌ **Gene list too large!** Your input contains {gene_set.size} genes, but the maximum allowed is 800 genes. Please reduce your gene list size.")
            raise typer.Exit(code=1)
        
        # Create unique results directory for this run with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        unique_folder_name = f"{gene_set_name}_{timestamp}"
        
        # Use output_dir if provided, otherwise use default results/ directory
        if output_dir is not None:
            gene_set_output_dir = output_dir / unique_folder_name
        else:
            gene_set_output_dir = ROOT / "results" / unique_folder_name
        
        gene_set_output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created unique output directory: {unique_folder_name}")
        
        if mode == "regular":
            run_regular_enrichment(
                gene_set, gene_set_libraries, background_gene_set,
                p_value_method, p_threshold, min_overlap,
                gene_set_output_dir
            )
        elif mode == "iterative":
            # Pass output_dir (may be None) to IterativeEnrichment for intermediate files
            # Pass gene_set_output_dir for final combined results
            run_iterative_enrichment(
                gene_set, gene_set_libraries, background_gene_set,
                p_value_method, p_threshold, min_overlap, min_term_size, max_term_size,
                max_iterations, output_dir, gene_set_output_dir, compute_benchmark=benchmark
            )


def run_regular_enrichment(
    gene_set: GeneSet,
    gene_set_libraries: List[GeneSetLibrary],
    background_gene_set: BackgroundGeneSet,
    p_value_method: str,
    p_threshold: float,
    min_overlap: int,
    output_dir: Path,
):
    """Run regular enrichment analysis."""
    logger.info("Running regular enrichment analysis")
    
    all_results = []
    
    for library in gene_set_libraries:
        logger.info(f"Processing library: {library.name}")
        
        try:
            enrich = Enrichment(
                gene_set,
                library,
                background_gene_set,
                min_term_size=library.library[0]["size"] if library.library else 10,
                max_term_size=library.library[-1]["size"] if library.library else 600,
                p_value_method_name=p_value_method,
            )
            
            # Filter results by p-value threshold and minimum overlap
            filtered_results = [
                result for result in enrich.results
                if (result.get("p-value", 1.0) <= p_threshold and
                    result.get("overlap_size", "").split("/")[0].isdigit() and 
                    int(result.get("overlap_size", "").split("/")[0]) >= min_overlap)
            ]
            
            # Add library name to each result for combined results
            for result in filtered_results:
                result["library"] = library.name
            
            # Update the enrichment object with filtered results to maintain compatibility
            enrich.results = filtered_results
            
            # Use the same formatting as Streamlit app
            if filtered_results:
                # Get the properly formatted DataFrame using the same method as Streamlit
                results_df = enrich.to_dataframe()
                output_file = output_dir / f"{library.name}_regular_results.tsv"
                results_df.to_csv(output_file, sep='\t', index=False)
                logger.info(f"Saved {len(filtered_results)} results to {output_file}")
                
                # Store for combined results (keep original format for compatibility)
                all_results.extend(filtered_results)
            
        except Exception as e:
            logger.error(f"Error processing library {library.name}: {e}")
            continue
    
    # Save combined results
    if all_results:
        # Create combined results in the same format as Streamlit
        combined_rows = []
        for result in all_results:
            import math
            # Get the library name from the result
            library_name = result.get("library", "")
            # If library name is not in the result, try to get it from the enrichment object
            if not library_name and hasattr(result, 'gene_set_library'):
                library_name = result.gene_set_library.name
            
            combined_rows.append({
                "Library": library_name,
                "Rank": result.get("rank", ""),
                "Term": result.get("term", "").replace("_", " ") if result.get("term") else "",
                "Description": result.get("description", ""),
                "Overlap size": result.get("overlap_size", ""),
                "p-value": result.get("p-value", ""),
                "-log(p-value)": -math.log10(result.get("p-value", 1)) if result.get("p-value", 1) > 0 else 0,
                "FDR": result.get("fdr", ""),
                "Genes": ", ".join(result.get("overlap", [])) if result.get("overlap") else "",
            })
        
        # Reorder columns to match Streamlit format
        column_order = ["Library", "Rank", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "FDR", "Genes"]
        combined_df = pd.DataFrame(combined_rows)
        if not combined_df.empty:
            combined_df = combined_df[column_order]
        
        combined_file = output_dir / "combined_regular_results.tsv"
        combined_df.to_csv(combined_file, sep='\t', index=False)
        logger.info(f"Saved combined results to {combined_file}")
        
        # Metadata is included in the text report files


def _get_library_display_name(library_name: str) -> str:
    """
    Return library display name. Since library.name is now always from alias.json,
    this function just returns it as-is for consistency.
    
    Args:
        library_name: Library name (should already be from alias.json "name" field)
        
    Returns:
        The library name (should already be the display name from alias.json)
    """
    # Since we now always use names from alias.json when creating GeneSetLibrary,
    # the library.name should already be the display name from alias.json
    return library_name


def _combine_dot_files(all_iter_results: dict) -> str:
    """
    Combine DOT files from multiple libraries into a single DOT file with colors.
    Uses the same coloring mechanism as the Streamlit app.
    
    Args:
        all_iter_results: Dictionary with library names (filenames) as keys and dict with 'results' and 'dot_content' as values
        
    Returns:
        Combined DOT content as string with colors applied
    """
    try:
        from ui.dot_utils import merge_iterative_dot
        
        # Map library names to display names for color matching
        # Since library.name is now always from alias.json, we can use it directly
        per_lib_dots = {}
        for lib_name, data in all_iter_results.items():
            # Library name is already from alias.json, use it directly
            per_lib_dots[lib_name] = data['dot_content']
        
        # Use the same merge function as Streamlit app
        combined_dot = merge_iterative_dot(per_lib_dots)
        return combined_dot
        
    except ImportError as e:
        logger.warning(f"Could not import merge_iterative_dot, using simple combination: {e}")
        # Fallback to simple combination without colors
        combined_dot = "graph iterative_enrichment {\n"
        combined_dot += "  graph [layout=neato];\n"
        combined_dot += "  node [shape=ellipse];\n"
        
        # Collect all nodes and edges
        all_nodes = set()
        all_edges = set()
        
        for lib_name, data in all_iter_results.items():
            dot_content = data['dot_content']
            
            # Parse the DOT content to extract nodes and edges
            lines = dot_content.split('\n')
            for line in lines:
                line = line.strip()
                if line.startswith('"') and '[' in line and ']' in line:
                    # This is a node definition
                    all_nodes.add(line)
                elif line.startswith('"') and ' -- ' in line:
                    # This is an edge definition
                    all_edges.add(line)
        
        # Add all nodes
        for node in sorted(all_nodes):
            combined_dot += f"  {node}\n"
        
        # Add all edges
        for edge in sorted(all_edges):
            combined_dot += f"  {edge}\n"
        
        # Close the graph
        combined_dot += "}\n"
        
        return combined_dot


def compute_benchmark_for_cli(
    iter_enrich_objects: Dict[str, IterativeEnrichment],
    gene_list_size: int,
    p_threshold: float,
    user_max_iterations: Optional[int],
    output_dir: Path,
    gene_set_name: str,
):
    """
    Compute statistical benchmarks for CLI (non-threaded version).
    
    Args:
        iter_enrich_objects: Dictionary of IterativeEnrichment objects by library name
        gene_list_size: Size of the input gene list
        p_threshold: P-value threshold used for enrichment
        user_max_iterations: User's max iterations (will be capped at 30)
        output_dir: Output directory for benchmark results
        gene_set_name: Name of the gene set
    """
    from parallel_null_distribution import (
        find_intersection_libraries,
        compute_null_distribution_parallel,
        normalize_library_name
    )
    from network_connectivity_benchmark import NetworkConnectivityAnalyzer
    from ui.benchmarking import generate_statistical_report_text, extract_benchmark_table_data
    
    # Find permutation data directory
    ROOT = Path(__file__).resolve().parent.parent
    parquet_dir = ROOT / "permutations" / "permutation_cluster_statistics_parquet"
    
    if not parquet_dir.exists():
        logger.warning(f"Permutation data directory not found: {parquet_dir}")
        logger.warning("Statistical benchmarking requires permutation data. Skipping benchmarks.")
        return
    
    # Get all library names from iter_enrich_objects
    # NOTE: All user-selected libraries are shown in enrichment results.
    # However, for statistical benchmarking, we only use libraries that match
    # between user selection and available permutation data in parquet files.
    user_selected_libraries = list(iter_enrich_objects.keys())
    
    # Find which user-selected libraries have matching permutation data
    # Statistics will only use libraries that match exactly between user selection and parquet files
    libraries_with_data, libraries_without_data = find_intersection_libraries(
        user_selected_libraries,
        parquet_dir,
        use_all_available=False  # Match user-selected libraries with available libraries
    )
    
    # Require minimum of 3 matching libraries for benchmarking
    # This ensures statistical robustness of the null distribution comparison
    MIN_MATCHING_LIBRARIES = 3
    if len(libraries_with_data) < MIN_MATCHING_LIBRARIES:
        logger.warning(
            f"Benchmarking requires at least {MIN_MATCHING_LIBRARIES} matching libraries. "
            f"Found {len(libraries_with_data)} matching libraries: {libraries_with_data}. "
            f"User selected: {user_selected_libraries}. "
            f"Skipping benchmarks."
        )
        return
    
    if not libraries_with_data:
        logger.warning("No libraries with permutation data available for benchmarking")
        logger.warning(f"Libraries in analysis: {user_selected_libraries}")
        logger.warning("Skipping benchmarks.")
        return
    
    logger.info(f"Found permutation data for {len(libraries_with_data)} libraries: {libraries_with_data}")
    if libraries_without_data:
        logger.warning(f"Libraries without permutation data (excluded from benchmarks): {libraries_without_data}")
    
    # Require p-value threshold <= 0.01 for benchmarking
    # This ensures statistical rigor and matches the permutation data generation
    if p_threshold > 0.01:
        logger.warning(
            f"Benchmarking requires p-value threshold <= 0.01. "
            f"Your threshold: {p_threshold}. "
            f"Skipping benchmarks."
        )
        return
    
    # Cap max iterations at 30 (permutation data maximum)
    if user_max_iterations is not None:
        user_max_iter = min(user_max_iterations, 30)
    else:
        user_max_iter = 30
    
    # Try to find merged permutation file
    merged_permutation_file = None
    possible_paths = [
        parquet_dir.parent / "merged_permutation_results.tsv",
        ROOT / "permutations" / "merged_permutation_results.tsv",
        ROOT / "results" / "permutation_results" / "merged_permutation_results.tsv",
        ROOT / "merged_permutation_results.tsv",
    ]
    for path in possible_paths:
        if path.exists():
            merged_permutation_file = path
            logger.info(f"Found merged permutation file: {merged_permutation_file}")
            break
    
    # Compute null distribution (non-threaded version)
    logger.info("Computing null distribution from permutation data...")
    null_dist_result = {
        'null_distribution': None,
        'status': 'running',
        'libraries_used': libraries_with_data,
        'error': None
    }
    import threading
    null_dist_lock = threading.Lock()
    
    # Call the parallel function directly (it will handle threading internally)
    compute_null_distribution_parallel(
        parquet_dir,
        gene_list_size,
        libraries_with_data,
        null_dist_result,
        null_dist_lock,
        p_threshold,
        user_max_iter,
        merged_permutation_file
    )
    
    # Check status
    final_status = null_dist_result.get('status', 'unknown')
    error_msg = null_dist_result.get('error', 'No error message')
    
    if final_status == 'error':
        logger.error(f"Failed to compute null distribution: {error_msg}")
        return
    elif final_status == 'unavailable':
        logger.warning(f"Statistical benchmarking unavailable: {error_msg}")
        return
    elif final_status != 'completed':
        logger.error(f"Unexpected status '{final_status}': {error_msg}")
        return
    
    null_distribution = null_dist_result.get('null_distribution')
    actual_size_used = null_dist_result.get('null_distribution_size', gene_list_size)
    
    if not null_distribution:
        logger.error("Null distribution is None")
        return
    
    # Verify the actual size is in the distribution
    actual_size_key = str(actual_size_used)
    if actual_size_key not in null_distribution:
        logger.warning(f"Size {actual_size_used} not found in null_distribution")
        size_keys = sorted([int(k) for k in null_distribution.keys()])
        if size_keys:
            actual_size_key = str(size_keys[0])
            actual_size_used = int(size_keys[0])
            logger.info(f"Using available size {actual_size_used} instead")
        else:
            logger.error("No sizes available in null distribution")
            return
    
    logger.info(f"Successfully computed null distribution for size {actual_size_used}")
    
    # Build combined network from ALL libraries to find clusters
    # (Clusters need all connections, but benchmarking will only use libraries with permutation data)
    combined_analyzer = NetworkConnectivityAnalyzer()
    
    # Determine iteration filter (cap at 30 to match permutation data)
    max_iter_filter = 30  # Always cap at 30 for fair comparison with permutation data
    if user_max_iter is not None:
        max_iter_filter = min(user_max_iter, 30)
        logger.info(f"Filtering iGEA results to iteration <= {max_iter_filter} (user requested {user_max_iter}, capped at 30)")
    else:
        # Check if user's results have more than 30 iterations
        max_user_iter = 0
        for iter_enrich in iter_enrich_objects.values():
            if iter_enrich.results:
                for record in iter_enrich.results:
                    iteration = record.get("Iteration", 1)
                    if isinstance(iteration, (int, float)):
                        max_user_iter = max(max_user_iter, int(iteration))
        if max_user_iter > 30:
            logger.info(f"User's results have {max_user_iter} iterations, trimming to 30 for fair comparison with permutation data")
    
    # Add results from ALL libraries to build the full network
    logger.info(f"Building network from all {len(iter_enrich_objects)} libraries to find clusters")
    for lib_name, iter_enrich in iter_enrich_objects.items():
        results = []
        for record in iter_enrich.results:
            # Filter by iteration count (always cap at 30)
            iteration = record.get("Iteration", 1)
            # Convert to int if needed
            if isinstance(iteration, str):
                try:
                    iteration = int(iteration)
                except (ValueError, TypeError):
                    iteration = 1
            elif not isinstance(iteration, (int, float)):
                iteration = 1
            else:
                iteration = int(iteration)
            
            # Filter: iteration must be <= max_iter_filter (which is always <= 30)
            if iteration > max_iter_filter:
                continue
            
            # Filter by p-value threshold
            # Try both 'p-value' and 'iteration p-value' field names
            p_value = record.get('p-value', record.get('iteration p-value', 1.0))
            if isinstance(p_value, str):
                try:
                    p_value = float(p_value)
                except (ValueError, TypeError):
                    p_value = 1.0
            if not isinstance(p_value, (int, float)):
                p_value = 1.0
            if p_value > p_threshold:
                continue
            
            # Get genes removed for this iteration (genes in the term)
            genes_removed = record.get("Genes removed for next iteration", [])
            if isinstance(genes_removed, str):
                genes_removed = [g.strip() for g in genes_removed.split(",") if g.strip()]
            elif not isinstance(genes_removed, list):
                # Fallback: try to get from 'Genes' field
                genes = record.get('Genes', [])
                if isinstance(genes, str):
                    genes_removed = [g.strip() for g in genes.split(",") if g.strip()]
                elif isinstance(genes, list):
                    genes_removed = genes
                else:
                    continue
            
            results.append({
                'Term': f"{lib_name}: {record.get('Term', '')}",
                'Iteration': iteration,
                'Library': lib_name,
                'Genes removed for next iteration': genes_removed,
            })
        
        if results:
            combined_analyzer.add_igea_results(results)
    
    # Create a mapping from user library names to parquet library names (for benchmarking)
    user_to_parquet_mapping = {}
    for user_lib in iter_enrich_objects.keys():
        normalized_user = normalize_library_name(user_lib)
        for parquet_lib in libraries_with_data:
            normalized_parquet = normalize_library_name(parquet_lib)
            if (normalized_user == normalized_parquet or 
                normalized_user in normalized_parquet.lower() or 
                normalized_parquet in normalized_user.lower()):
                user_to_parquet_mapping[user_lib] = parquet_lib
                logger.info(f"Matched user library '{user_lib}' to parquet library '{parquet_lib}'")
                break
    
    if not user_to_parquet_mapping:
        logger.warning("No matching libraries found between analysis and permutation data")
        logger.warning(f"User libraries: {list(iter_enrich_objects.keys())}")
        logger.warning(f"Parquet libraries: {libraries_with_data}")
        return
    
    logger.info(f"Network built from all libraries. Benchmarking will use {len(user_to_parquet_mapping)} libraries with permutation data")
    
    # Get clusters
    clusters = combined_analyzer.get_clusters()
    
    if not clusters:
        logger.warning("No clusters found in the network. Benchmarking requires at least one cluster.")
        # Generate a text report explaining the situation
        report_text = f"""Statistical Benchmarking Report
{'=' * 80}

Gene Set: {gene_set_name}
Gene Set Size: {gene_list_size}
Actual Size Used for Null Distribution: {actual_size_used}
P-value Threshold: {p_threshold}
Max Iterations: {user_max_iter if user_max_iter else 'Unlimited'}

Libraries with Permutation Data: {len(libraries_with_data)}
{', '.join(libraries_with_data) if libraries_with_data else 'None'}

Libraries without Permutation Data: {len(libraries_without_data)}
{', '.join(libraries_without_data) if libraries_without_data else 'None'}

{'=' * 80}

STATUS: No Clusters Found

Benchmarking requires at least one cluster in the network to compare against 
the null distribution. A cluster is formed when enriched terms from different 
libraries share common genes, creating a connected subgraph in the network.

In this analysis, no clusters were detected, which means:
- The enriched terms from different libraries do not share enough common genes
- The network consists of isolated nodes (terms and genes) without connections
- The iterative enrichment process found terms, but they don't form connected clusters

This can occur when:
- The gene list has diverse functions that don't overlap across libraries
- The p-value threshold is too strict, resulting in few enriched terms
- The libraries used have minimal overlap in their gene sets

Null Distribution Summary:
The null distribution was successfully computed from {len(null_distribution)} gene list sizes.
However, without clusters to benchmark, statistical comparison cannot be performed.

{'=' * 80}
Report generated: {datetime.now().isoformat()}
"""
        
        report_file = output_dir / f"{gene_set_name}_statistical_report.txt"
        with open(report_file, 'w') as f:
            f.write(report_text)
        logger.info(f"Saved statistical report (no clusters) to {report_file}")
        
        # Create an empty table file to indicate benchmarking was attempted
        table_file = output_dir / "statistical_benchmarks_table.tsv"
        table_df = pd.DataFrame([{
            "Status": "No clusters found",
            "Message": "Benchmarking requires at least one cluster in the network"
        }])
        table_df.to_csv(table_file, sep='\t', index=False)
        logger.info(f"Saved benchmark table (no clusters) to {table_file}")
        
        logger.warning("Benchmarking files generated but no clusters were found to benchmark")
        return
    
    logger.info(f"Found {len(clusters)} clusters in the network")
    
    # Benchmark only the largest cluster (null distribution is built from largest clusters only)
    # Clusters are already sorted by size (largest first)
    largest_cluster = clusters[0]
    
    logger.info(f"Benchmarking largest cluster: {largest_cluster['size']} nodes "
               f"({largest_cluster['n_genes']} genes, {largest_cluster['n_terms']} terms)")
    logger.info(f"Note: Only the largest cluster is benchmarked against null distribution "
               f"(which is built from largest clusters in each permutation)")
    
    # Benchmark the largest cluster
    # But only if it contains at least one library with permutation data
    cluster_benchmarks = []
    
    # Check if the largest cluster contains any libraries with permutation data
    # Get libraries from the analyzer's term_to_library mapping (most reliable)
    cluster_libraries = set()
    for term in largest_cluster.get('terms', []):
        # Get library from the analyzer's term_to_library mapping
        if term in combined_analyzer.term_to_library:
            lib_name = combined_analyzer.term_to_library[term]
            cluster_libraries.add(lib_name)
        else:
            # Fallback: Extract library name from term (format: "library_name: term_name")
            # Try to match against known library names first (to handle full names with colons)
            matched = False
            for lib_name in iter_enrich_objects.keys():
                if term.startswith(f"{lib_name}:"):
                    cluster_libraries.add(lib_name)
                    matched = True
                    break
            if not matched and ':' in term:
                # If no match, just take the first part before colon (less reliable)
                lib_name = term.split(':', 1)[0]
                cluster_libraries.add(lib_name)
    
    # Check if any cluster library matches libraries with permutation data
    cluster_has_benchmarkable_libs = False
    for cluster_lib in cluster_libraries:
        normalized_cluster = normalize_library_name(cluster_lib)
        for parquet_lib in libraries_with_data:
            normalized_parquet = normalize_library_name(parquet_lib)
            if (normalized_cluster == normalized_parquet or 
                normalized_cluster in normalized_parquet.lower() or 
                normalized_parquet in normalized_cluster.lower()):
                cluster_has_benchmarkable_libs = True
                break
        if cluster_has_benchmarkable_libs:
            break
    
    if cluster_has_benchmarkable_libs:
        from network_connectivity_benchmark import benchmark_cluster
        benchmark = benchmark_cluster(
            largest_cluster,
            null_distribution,
            actual_size_used,
            use_interpolation=False
        )
        if benchmark:
            # Convert sets to lists for JSON serialization
            cluster_for_json = {
                'cluster_number': largest_cluster.get('cluster_number'),
                'size': largest_cluster.get('size'),
                'n_genes': largest_cluster.get('n_genes'),
                'n_terms': largest_cluster.get('n_terms'),
                'n_edges': largest_cluster.get('n_edges'),
                'density': largest_cluster.get('density'),
                'avg_edges_per_gene': largest_cluster.get('avg_edges_per_gene'),
                'n_libraries': largest_cluster.get('n_libraries'),
                'genes': sorted(list(largest_cluster.get('genes', []))),
                'terms': sorted(list(largest_cluster.get('terms', []))),
            }
            cluster_benchmarks.append({
                'cluster': cluster_for_json,
                'benchmark': benchmark
            })
        else:
            logger.warning(f"benchmark_cluster returned empty result for largest cluster "
                         f"with {largest_cluster.get('n_genes', 0)} genes, {largest_cluster.get('n_terms', 0)} terms")
    else:
        logger.warning(f"Largest cluster does not contain any libraries with permutation data. "
                      f"Cluster libraries: {cluster_libraries}, "
                      f"Libraries with data: {libraries_with_data}")
    
    logger.info(f"Computed benchmarks for {len(cluster_benchmarks)} clusters")
    
    # Save benchmark results
    # Generate and save text report
    report_text = generate_statistical_report_text(
        cluster_benchmarks,
        gene_set_name,
        libraries_with_data,
        libraries_without_data,
        combined_analyzer
    )
    
    report_file = output_dir / f"{gene_set_name}_statistical_report.txt"
    with open(report_file, 'w') as f:
        f.write(report_text)
    logger.info(f"Saved statistical report to {report_file}")
    
    # 3. Save benchmark table as TSV
    table_data = extract_benchmark_table_data(cluster_benchmarks)
    if table_data:
        table_df = pd.DataFrame(table_data)
        table_file = output_dir / "statistical_benchmarks_table.tsv"
        table_df.to_csv(table_file, sep='\t', index=False)
        logger.info(f"Saved benchmark table to {table_file}")
    
    logger.info("✅ Statistical benchmarking completed successfully!")


def run_iterative_enrichment(
    gene_set: GeneSet,
    gene_set_libraries: List[GeneSetLibrary],
    background_gene_set: BackgroundGeneSet,
    p_value_method: str,
    p_threshold: float,
    min_overlap: int,
    min_term_size: int,
    max_term_size: int,
    max_iterations: int,
    output_dir: Optional[Path],  # For IterativeEnrichment intermediate files (None = use default results/)
    final_output_dir: Path,  # For final combined results
    compute_benchmark: bool = False,
):
    """Run iterative enrichment analysis."""
    logger.info("Running iterative enrichment analysis")
    
    # Generate a shared run ID for all libraries
    shared_run_id = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]
    logger.info(f"Generated shared run ID: {shared_run_id}")
    
    all_iter_results = {}
    all_iter_enrich_objects = {}  # Store IterativeEnrichment objects for benchmarking
    
    for library in gene_set_libraries:
        logger.info(f"Processing library: {library.name}")
        
        try:
            # Create progress callback for logging
            def progress_callback(message: str):
                logger.info(f"Library {library.name}: {message}")
            
            # Pass output_dir to IterativeEnrichment
            # If None, IterativeEnrichment will use its default results/ directory
            # If provided, it will use that directory for intermediate files
            it = IterativeEnrichment(
                gene_set=gene_set,
                gene_set_library=library,
                background_gene_set=background_gene_set,
                min_term_size=min_term_size,
                max_term_size=max_term_size,
                p_value_method_name=p_value_method,
                p_threshold=p_threshold,
                max_iterations=None if max_iterations == 0 else max_iterations,
                min_overlap=min_overlap,
                progress_callback=progress_callback,
                run_id=shared_run_id,
                output_dir=output_dir,  # None = use default results/ directory
            )
            
            # Save individual library results
            if it.results:
                # Use the same formatting as Streamlit app
                results_df = it.to_dataframe()
                
                output_file = final_output_dir / f"{library.name}_iterative_results.tsv"
                results_df.to_csv(output_file, sep='\t', index=False)
                logger.info(f"Saved {len(it.results)} iterations to {output_file}")
                
                # Store DOT content for combined network
                dot_content = it.to_dot()
                all_iter_results[library.name] = {
                    'results': it.results,
                    'dot_content': dot_content
                }
                
                # Store IterativeEnrichment object for benchmarking if needed
                if compute_benchmark:
                    all_iter_enrich_objects[library.name] = it
            
        except Exception as e:
            logger.error(f"Error processing library {library.name}: {e}")
            continue
    
    # Save combined results
    if all_iter_results:
        # Create combined TSV in the same format as Streamlit
        combined_rows = []
        for lib_name, data in all_iter_results.items():
            results = data['results']
            for result in results:
                import math
                combined_rows.append({
                    "Library": lib_name,
                    "Iteration": result.get("Iteration", ""),
                    "Term": result.get("Term", "").replace("_", " ") if result.get("Term") else "",
                    "Description": result.get("Description", ""),
                    "Overlap size": result.get("Overlap size", ""),
                    "p-value": result.get("p-value", ""),
                    "-log(p-value)": -math.log10(result.get("p-value", 1)) if result.get("p-value", 1) > 0 else 0,
                    "Genes": ", ".join(result.get('Genes', [])) if result.get('Genes') else "",
                })
        
        if combined_rows:
            # Reorder columns to match Streamlit format
            column_order = ["Library", "Iteration", "Term", "Description", "Overlap size", "p-value", "-log(p-value)", "Genes"]
            combined_df = pd.DataFrame(combined_rows)
            if not combined_df.empty:
                combined_df = combined_df[column_order]
            
            combined_file = final_output_dir / "combined_iterative_results.tsv"
            combined_df.to_csv(combined_file, sep='\t', index=False)
            logger.info(f"Saved combined iterative results to {combined_file}")
        
        # Generate combined DOT file
        try:
            combined_dot_content = _combine_dot_files(all_iter_results)
            combined_dot_file = final_output_dir / "combined_network.dot"
            with open(combined_dot_file, 'w') as f:
                f.write(combined_dot_content)
            logger.info(f"Saved combined network DOT to {combined_dot_file}")
            
            # Generate and save AI analysis prompt for combined network
            try:
                from ui.rendering import generate_ai_analysis_prompt
                
                # Enhanced prompt format
                enhanced_prompt = generate_ai_analysis_prompt(combined_dot_content)
                enhanced_file = final_output_dir / "combined_ai_analysis_prompt_enhanced.txt"
                with open(enhanced_file, 'w') as f:
                    f.write(enhanced_prompt)
                logger.info(f"Saved combined enhanced AI prompt to {enhanced_file}")
                
            except ImportError as e:
                logger.warning(f"Could not import AI analysis functions: {e}")
            except Exception as e:
                logger.warning(f"Error generating AI prompts: {e}")
                
        except Exception as e:
            logger.warning(f"Error generating combined DOT file: {e}")
        
        # Metadata is included in the text report files
        
        # Compute statistical benchmarks if requested
        if compute_benchmark and all_iter_enrich_objects:
            logger.info("Computing statistical benchmarks...")
            try:
                compute_benchmark_for_cli(
                    all_iter_enrich_objects,
                    gene_set.size,
                    p_threshold,
                    max_iterations if max_iterations > 0 else None,
                    final_output_dir,
                    gene_set.name
                )
            except Exception as e:
                logger.error(f"Error computing benchmarks: {e}")
                logger.warning("Benchmarking failed, but enrichment analysis completed successfully")


@app.command(help="Run gene set enrichment analysis from the command line")
def main(
    gene_sets: List[Path] = typer.Option(
        None,
        "--genelist",
        "-g",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Paths to gene set files.",
    ),
    gene_format: str = typer.Option(
        "symbols",
        "--gene-format",
        help="Gene input format: 'symbols' (Gene Symbols) or 'entrez_ids' (Entrez IDs)",
    ),
    background: Optional[Path] = typer.Option(
        None,
        "--background",
        "-b",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Path to the background gene set file. Default: all_genes.txt",
        show_default=False,
    ),
    libraries: List[Path] = typer.Option(
        None,
        "--libraries",
        "-l",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Paths to gene set library files. Default: active libraries from alias.json.",
        show_default=False,
    ),
    mode: str = typer.Option(
        "regular",
        "--mode",
        "-m",
        help="Analysis mode: 'regular' or 'iterative'",
    ),
    p_value_method: str = typer.Option(
        "Fisher's Exact Test",
        "--method",
        help="P-value calculation method: 'fisher' (Fisher's Exact Test), 'hga' (Hypergeometric Test), or 'chi' (Chi-squared Test)",
    ),
    p_threshold: float = typer.Option(
        0.01,
        "--p-threshold",
        "-p",
        help="Raw p-value threshold for including terms (not FDR-corrected)",
    ),
    min_overlap: int = typer.Option(
        3,
        "--min-overlap",
        help="Minimum overlap size required for terms",
    ),
    min_term_size: int = typer.Option(
        10,
        "--min-term-size",
        help="Minimum term size",
    ),
    max_term_size: int = typer.Option(
        600,
        "--max-term-size",
        help="Maximum term size",
    ),
    max_iterations: int = typer.Option(
        10,
        "--max-iterations",
        help="Maximum iterations for iterative mode (0 = no limit)",
    ),
    output_dir: Optional[Path] = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Output directory for results",
    ),
    benchmark: bool = typer.Option(
        False,
        "--benchmark",
        help="Compute statistical benchmarks against permutation data (iterative mode only)",
    ),
):
    """
    Run gene set enrichment analysis from the command line.
    
    This CLI provides the same core functionality as the Streamlit app,
    supporting both regular and iterative enrichment analysis modes.
    """
    
    # Convert method shortcuts to full names
    method_mapping = {
        "fisher": "Fisher's Exact Test",
        "hga": "Hypergeometric Test", 
        "chi": "Chi-squared Test"
    }
    
    if p_value_method.lower() in method_mapping:
        p_value_method = method_mapping[p_value_method.lower()]
    
    # Validate gene format
    if gene_format not in ["symbols", "entrez_ids"]:
        typer.echo("Error: Gene format must be 'symbols' or 'entrez_ids'", err=True)
        raise typer.Exit(code=1)
    
    # Validate mode
    if mode not in ["regular", "iterative"]:
        typer.echo("Error: Mode must be 'regular' or 'iterative'", err=True)
        raise typer.Exit(code=1)
    
    # Default values handling
    if not gene_sets:
        gene_sets = list((ROOT / "data/gene_lists").glob("*.txt"))
        if not gene_sets:
            typer.echo("Error: No gene set files found in data/gene_lists/", err=True)
            raise typer.Exit(code=1)

    if not background:
        # Try to use all_genes.txt as default
        all_genes_file = ROOT / "data/backgrounds/all_genes.txt"
        if all_genes_file.exists():
            background = all_genes_file
            typer.echo(f"🎯 Using default background: all_genes.txt")
        else:
            # Fallback to first available background file
            background_files = list((ROOT / "data/backgrounds").glob("*.txt"))
            background = background_files[0] if background_files else None
            if not background:
                typer.echo("Error: No background files found in data/backgrounds/", err=True)
                raise typer.Exit(code=1)
            typer.echo(f"⚠️  all_genes.txt not found, using: {background.name}")

    if not libraries:
        # Load active libraries from alias file
        alias_file = ROOT / "data/libraries/alias.json"
        if not alias_file.exists():
            typer.echo("Error: alias.json file not found in data/libraries/", err=True)
            raise typer.Exit(code=1)
        
        try:
            import json
            with open(alias_file, 'r') as f:
                aliases = json.load(f)
            
            # Get active libraries
            active_libraries = []
            for alias in aliases:
                if alias.get("active", False):
                    lib_file = ROOT / "data/libraries" / alias["file"]
                    if lib_file.exists():
                        active_libraries.append(lib_file)
                    else:
                        typer.echo(f"Warning: Active library file not found: {alias['file']}", err=True)
            
            if not active_libraries:
                typer.echo("Error: No active library files found in alias.json", err=True)
                raise typer.Exit(code=1)
            
            libraries = active_libraries
            typer.echo(f"📚 Using active libraries from alias.json: {len(libraries)} files")
            for lib in libraries:
                typer.echo(f"   - {lib.name}")
            typer.echo("")
            
        except Exception as e:
            typer.echo(f"Error reading alias.json: {e}", err=True)
            raise typer.Exit(code=1)

    # Ensure that the resulting variables are not empty
    if not gene_sets or not libraries:
        typer.echo("Error: Gene sets and libraries cannot be empty.", err=True)
        raise typer.Exit(code=1)

    # Display analysis parameters
    typer.echo(f"Analysis Mode: {mode}")
    typer.echo(f"Gene Sets: {len(gene_sets)} files")
    typer.echo(f"Background: {background}")
    typer.echo(f"Libraries: {len(libraries)} files")
    typer.echo(f"P-value Method: {p_value_method}")
    typer.echo(f"P-value Threshold: {p_threshold}")
    typer.echo(f"Min Overlap: {min_overlap}")
    typer.echo(f"Term Size Range: {min_term_size}-{max_term_size}")
    if mode == "iterative":
        typer.echo(f"Max Iterations: {max_iterations}")
        typer.echo(f"Statistical Benchmarking: {'Enabled' if benchmark else 'Disabled'}")
    # Determine output directory
    # If not provided, use None to let IterativeEnrichment use default results/ directory
    # If provided, use it for both final results and IterativeEnrichment's internal directory
    final_output_dir = output_dir
    if output_dir is not None:
        typer.echo(f"Output Directory: {output_dir}")
    else:
        typer.echo(f"Output Directory: {ROOT / 'results'} (default)")
    typer.echo("")
    
    # Run the enrichment analysis
    try:
        run_enrichment(
            gene_sets=gene_sets,
            background=background,
            libraries=libraries,
            gene_format=gene_format,
            p_value_method=p_value_method,
            mode=mode,
            p_threshold=p_threshold,
            min_overlap=min_overlap,
            min_term_size=min_term_size,
            max_term_size=max_term_size,
            max_iterations=max_iterations,
            output_dir=final_output_dir,  # None = use default results/ directory
            benchmark=benchmark,
        )
        if final_output_dir is not None:
            typer.echo(f"✅ Analysis completed successfully! Results saved to {final_output_dir}")
        else:
            typer.echo(f"✅ Analysis completed successfully! Results saved to {ROOT / 'results'}")
    except Exception as e:
        typer.echo(f"❌ Error during analysis: {e}", err=True)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
