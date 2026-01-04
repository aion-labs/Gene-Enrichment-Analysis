#!/usr/bin/env python3
"""
Benchmark HIV gene list network connectivity against null distribution.

This script:
1. Loads HIV gene list
2. Runs iGEA for 11 libraries (all libraries with permutation data)
3. Computes network connectivity metrics
4. Benchmarks against null distribution from permutations
5. Displays results
"""

import sys
import logging
from pathlib import Path
from typing import Dict, List
import json
import pandas as pd

# Add code directory to path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from iter_enrichment import IterativeEnrichment
from gene_converter import GeneConverter
from network_connectivity_benchmark import (
    NetworkConnectivityAnalyzer,
    benchmark_real_results
)
from parallel_null_distribution import (
    find_intersection_libraries,
    compute_null_distribution_parallel,
    normalize_library_name
)
import threading

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Paths
DATA_DIR = PROJECT_ROOT / "data"
LIBRARIES_DIR = DATA_DIR / "libraries"
BACKGROUNDS_DIR = DATA_DIR / "backgrounds"
GENE_LISTS_DIR = DATA_DIR / "gene_lists"
NULL_DIST_FILE = PROJECT_ROOT / "results" / "connectivity_null_distribution.json"
PARQUET_DIR = PROJECT_ROOT / "permutations" / "permutation_cluster_statistics_parquet"

# Libraries to use (matching permutation data - all 11 libraries)
# Note: Library names should match parquet library names for proper matching
LIBRARIES = {
    "Reactome": "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
    "KEGG": "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",
    "GO BP": "c5.go.bp.v2025.1.Hs.symbols.gmt",
    "GO MF": "c5.go.mf.v2025.1.Hs.symbols.gmt",
    "GO CC": "c5.go.cc.v2025.1.Hs.symbols.gmt",
    "BioCarta": "c2.cp.biocarta.v2025.1.Hs.symbols.gmt",
    "Canonical pathways": "c2.cp.v2025.1.Hs.symbols.gmt",
    "KEGG Medicus": "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt",
    "Pathway Interaction Database": "c2.cp.pid.v2025.1.Hs.symbols.gmt",
    "WikiPathways": "c2.cp.wikipathways.v2025.1.Hs.symbols.gmt",
    "Hallmark": "h.all.v2025.1.Hs.symbols.gmt",
}

# iGEA parameters
PARAMS = {
    "p_threshold": 0.01,  # Using 0.01 for benchmarking
    "min_overlap": 3,
    "min_term_size": 10,
    "max_term_size": 600,
    "max_iterations": 10,
    "p_value_method": "Fisher's Exact Test",
}


def load_hiv_gene_list() -> List[str]:
    """Load HIV gene list and convert Entrez IDs to gene symbols."""
    logger.info("Loading HIV gene list...")
    
    hiv_file = GENE_LISTS_DIR / "HIV.InputGeneList.txt"
    if not hiv_file.exists():
        raise FileNotFoundError(f"HIV gene list not found: {hiv_file}")
    
    # Read Entrez IDs
    with open(hiv_file, 'r') as f:
        entrez_ids = [line.strip() for line in f if line.strip()]
    
    logger.info(f"Loaded {len(entrez_ids)} Entrez IDs from HIV gene list")
    
    # Convert to gene symbols
    converter = GeneConverter()
    gene_symbols = []
    failed = []
    
    for entrez_id in entrez_ids:
        symbol = converter.get_symbol(entrez_id)
        if symbol:
            gene_symbols.append(symbol)
        else:
            failed.append(entrez_id)
    
    if failed:
        logger.warning(f"Failed to convert {len(failed)} Entrez IDs to symbols")
    
    logger.info(f"Converted to {len(gene_symbols)} gene symbols")
    return gene_symbols


def load_libraries() -> Dict[str, GeneSetLibrary]:
    """Load all 11 libraries for which we have permutation data."""
    logger.info("Loading gene set libraries...")
    libraries = {}
    
    for lib_name, filename in LIBRARIES.items():
        lib_path = LIBRARIES_DIR / filename
        if not lib_path.exists():
            logger.warning(f"Library file not found: {lib_path}, skipping {lib_name}")
            continue
        
        logger.info(f"Loading library: {lib_name}")
        library = GeneSetLibrary(str(lib_path), name=lib_name)
        
        # Filter terms by size
        filtered_terms = [
            t for t in library.library
            if PARAMS["min_term_size"] <= t["size"] <= PARAMS["max_term_size"]
        ]
        library.library = filtered_terms
        library.num_terms = len(filtered_terms)
        library.unique_genes = library.compute_unique_genes()
        library.size = len(library.unique_genes)
        
        libraries[lib_name] = library
        logger.info(f"  {lib_name}: {library.num_terms} terms, {library.size} unique genes")
    
    return libraries


def load_background() -> BackgroundGeneSet:
    """Load the background gene set."""
    bg_path = BACKGROUNDS_DIR / "all_genes.txt"
    if not bg_path.exists():
        raise FileNotFoundError(f"Background file not found: {bg_path}")
    
    logger.info(f"Loading background: {bg_path}")
    bg = BackgroundGeneSet(str(bg_path), name="all_genes", input_format="symbols", skip_validation=True)
    logger.info(f"Loaded background: {bg.size} genes")
    return bg


def run_igea_for_library(
    gene_set: GeneSet,
    library: GeneSetLibrary,
    background: BackgroundGeneSet
) -> IterativeEnrichment:
    """Run iGEA for a single library."""
    logger.info(f"Running iGEA for {library.name}...")
    
    iter_enrich = IterativeEnrichment(
        gene_set=gene_set,
        gene_set_library=library,
        background_gene_set=background,
        min_term_size=PARAMS["min_term_size"],
        max_term_size=PARAMS["max_term_size"],
        p_value_method_name=PARAMS["p_value_method"],
        p_threshold=PARAMS["p_threshold"],
        max_iterations=PARAMS["max_iterations"],
        min_overlap=PARAMS["min_overlap"],
        progress_callback=None,
        use_multiprocessing=True,
    )
    
    logger.info(f"  {library.name}: {len(iter_enrich.results)} iterations")
    return iter_enrich


# Removed: compute_connectivity_for_library - not needed since we only analyze combined network


def print_benchmark_results(benchmark: Dict, library_name: str):
    """Print benchmark results in a readable format."""
    print("\n" + "=" * 80)
    print(f"Library: {library_name}")
    print("=" * 80)
    
    if not benchmark or 'comparison' not in benchmark:
        print("No benchmark results available")
        return
    
    real_metrics = benchmark.get('real_metrics', {})
    comparison = benchmark.get('comparison', {})
    
    print(f"\nNetwork Summary:")
    print(f"  Genes: {real_metrics.get('n_genes', 0)}")
    print(f"  Terms: {real_metrics.get('n_terms', 0)}")
    print(f"  Edges: {real_metrics.get('n_edges', 0)}")
    
    # Show size information
    if 'size_difference' in benchmark:
        print(f"\nSize Matching:")
        print(f"  Real gene list size: {benchmark.get('gene_list_size', 'N/A')}")
        size_diff = benchmark.get('size_difference', 0)
        nearest_size = benchmark.get('nearest_null_size', 'N/A')
        if size_diff == 0:
            print(f"  ✓ Exact match with null distribution size: {nearest_size}")
        elif size_diff <= 50 and nearest_size > benchmark.get('gene_list_size', 0):
            # Likely rounded up to next increment of 50
            print(f"  ✓ Using null distribution size: {nearest_size} (rounded from {benchmark.get('gene_list_size', 'N/A')}, diff: {size_diff} genes)")
        else:
            print(f"  Interpolated between sizes: {nearest_size} (diff: {size_diff} genes, {benchmark.get('size_difference_pct', 0):.1f}%)")
            if size_diff > 25:
                print(f"  ⚠ Warning: Large size difference - results may be less accurate")
    
    print(f"\n{'Metric':<35} {'Real Value':<15} {'Null Mean':<15} {'Z-score':<12} {'Percentile':<12} {'Better?':<10}")
    print("-" * 100)
    
    metrics_order = [
        'avg_connections_per_gene',
        'network_density',
        'n_connected_components',
        'largest_component_size',
        'gene_centrality_max',
        'clustering_coefficient',
    ]
    
    for metric in metrics_order:
        if metric not in comparison:
            continue
        
        comp = comparison[metric]
        real_val = comp['real_value']
        null_mean = comp['null_mean']
        z_score = comp['z_score']
        percentile = comp['percentile']
        is_better = "✓ YES" if comp['is_better'] else "✗ NO"
        
        # Format based on metric type
        if metric == 'n_connected_components':
            # Lower is better for components
            is_better = "✓ YES" if real_val < null_mean else "✗ NO"
        
        print(f"{metric:<35} {real_val:<15.4f} {null_mean:<15.4f} {z_score:<12.2f} {percentile:<12.1f}% {is_better:<10}")
    
    # Summary
    better_count = sum(1 for comp in comparison.values() if comp['is_better'])
    total_count = len(comparison)
    
    # Special handling for n_connected_components (lower is better)
    if 'n_connected_components' in comparison:
        comp = comparison['n_connected_components']
        if comp['real_value'] < comp['null_mean']:
            better_count += 1
        total_count += 1
    
    print("\n" + "-" * 100)
    print(f"Summary: {better_count}/{total_count} metrics are better than null")
    
    if better_count == total_count:
        print("✓ Network connectivity is SIGNIFICANTLY BETTER than random!")
    elif better_count >= total_count * 0.7:
        print("⚠ Network connectivity is MODERATELY BETTER than random")
    else:
        print("✗ Network connectivity is NOT significantly better than random")


def main():
    """Main function."""
    print("\n" + "=" * 80)
    print("HIV Gene List - Network Connectivity Benchmarking")
    print("=" * 80)
    
    # 1. Load HIV gene list
    print("\n[1/5] Loading HIV gene list...")
    hiv_genes = load_hiv_gene_list()
    print(f"✓ Loaded {len(hiv_genes)} genes")
    
    # 2. Load background
    print("\n[2/5] Loading background gene set...")
    background = load_background()
    print(f"✓ Loaded background: {background.size} genes")
    
    # 3. Create gene set
    print("\n[3/5] Creating gene set...")
    gene_set = GeneSet(
        gene_list=hiv_genes,
        validation_set=background.genes,
        hgcn=False,
        format=False
    )
    print(f"✓ Gene set created: {gene_set.size} genes")
    
    # 4. Load libraries
    print("\n[4/6] Loading gene set libraries...")
    libraries = load_libraries()
    print(f"✓ Loaded {len(libraries)} libraries")
    
    # 5. Check which libraries have permutation data
    print("\n[5/6] Checking libraries with permutation data...")
    user_selected_libraries = list(libraries.keys())
    libraries_with_data, libraries_without_data = find_intersection_libraries(
        user_selected_libraries,
        PARQUET_DIR,
        use_all_available=False  # Match user-selected libraries with available libraries
    )
    
    # Require minimum of 3 matching libraries for benchmarking
    MIN_MATCHING_LIBRARIES = 3
    if len(libraries_with_data) < MIN_MATCHING_LIBRARIES:
        raise ValueError(
            f"Benchmarking requires at least {MIN_MATCHING_LIBRARIES} matching libraries. "
            f"Found {len(libraries_with_data)} matching libraries: {libraries_with_data}. "
            f"User selected: {user_selected_libraries}. "
            f"Please ensure at least {MIN_MATCHING_LIBRARIES} libraries match between your selection and parquet files."
        )
    
    if not libraries_with_data:
        raise ValueError(
            f"None of the selected libraries have permutation data available. "
            f"Selected: {user_selected_libraries}. "
            f"Please ensure Parquet files exist in {PARQUET_DIR}"
        )
    
    print(f"✓ Libraries with permutation data: {len(libraries_with_data)}")
    print(f"  {', '.join(libraries_with_data)}")
    if libraries_without_data:
        print(f"⚠ Libraries without permutation data (included in enrichment, excluded from statistics):")
        print(f"  {', '.join(libraries_without_data)}")
    
    # 6. Start null distribution computation in parallel thread
    print("\n[6/6] Starting parallel null distribution computation...")
    
    # Check if p-value threshold is valid for benchmarking
    user_p_threshold = PARAMS["p_threshold"]
    if user_p_threshold > 0.01:
        print(f"⚠ WARNING: Benchmarking requires p-value threshold <= 0.01")
        print(f"   Your p-value threshold: {user_p_threshold}")
        print(f"   Statistical benchmarking will be skipped.")
        null_dist_result = {
            'null_distribution': None,
            'status': 'unavailable',
            'libraries_used': [],
            'error': f"P-value threshold {user_p_threshold} > 0.01"
        }
    else:
        null_dist_result = {
            'null_distribution': None,
            'status': 'running',
            'libraries_used': libraries_with_data,
            'error': None
        }
        null_dist_lock = threading.Lock()
        
        # Try to find merged permutation file for p-value filtering
        merged_permutation_file = None
        possible_paths = [
            PROJECT_ROOT / "permutations" / "merged_permutation_results.tsv",
            PROJECT_ROOT / "results" / "permutation_results" / "merged_permutation_results.tsv",
            PROJECT_ROOT / "archive" / "permutation_analysis" / "results" / "permutation_results" / "merged_permutation_results.tsv",
            PROJECT_ROOT / "merged_permutation_results.tsv",
        ]
        for path in possible_paths:
            if path.exists():
                merged_permutation_file = path
                logger.info(f"Found merged permutation file: {merged_permutation_file}")
                break
        
        null_dist_thread = threading.Thread(
            target=compute_null_distribution_parallel,
            args=(
                PARQUET_DIR,
                gene_set.size,
                libraries_with_data,
                null_dist_result,
                null_dist_lock,
                user_p_threshold,  # Pass user's p-value threshold
                PARAMS["max_iterations"],  # Pass max iterations (capped at 30)
                merged_permutation_file  # Pass merged permutation file path
            ),
            daemon=True
        )
        null_dist_thread.start()
        print(f"✓ Null distribution computation started in parallel thread")
    
    # 7. Run iGEA for each library (only to combine later, no individual analysis)
    print("\n" + "=" * 80)
    print("Running iGEA for each library...")
    print("=" * 80)
    print(f"Note: All {len(libraries)} libraries will be included in enrichment results.")
    if libraries_without_data:
        print(f"      However, statistics will only be computed for {len(libraries_with_data)} libraries.")
    print()
    
    # Check if we can reuse existing results
    import os
    SKIP_IGEA = os.environ.get('SKIP_IGEA', 'false').lower() == 'true'
    
    all_results = {}
    
    if SKIP_IGEA:
        print("⚠ SKIP_IGEA environment variable set - attempting to load existing results...")
        # Try to load from existing benchmark JSON if it exists
        benchmark_file = PROJECT_ROOT / "results" / "hiv_connectivity_benchmark" / "hiv_connectivity_benchmark.json"
        if benchmark_file.exists():
            try:
                with open(benchmark_file, 'r') as f:
                    existing_data = json.load(f)
                print(f"✓ Found existing benchmark file: {benchmark_file}")
                print("  Note: This mode requires existing iGEA results. If results are missing, set SKIP_IGEA=false")
                # We still need to run iGEA to get the results, but we can skip if results exist
                # For now, we'll still run it but this gives us a hook for future optimization
            except Exception as e:
                logger.warning(f"Could not load existing benchmark: {e}")
                SKIP_IGEA = False
    
    for lib_name, library in libraries.items():
        try:
            iter_enrich = run_igea_for_library(gene_set, library, background)
            all_results[lib_name] = iter_enrich
            print(f"✓ {lib_name}: {len(iter_enrich.results)} iterations")
        except Exception as e:
            logger.error(f"Error processing {lib_name}: {e}", exc_info=True)
            all_results[lib_name] = None
    
    # Wait for null distribution computation to complete (if it was started)
    if null_dist_result['status'] == 'running':
        print("\n" + "=" * 80)
        print("Waiting for null distribution computation...")
        print("=" * 80)
        null_dist_thread.join(timeout=60)  # Wait up to 60 seconds (longer for raw data processing)
        
        # Re-check status after thread completes
        final_status = null_dist_result.get('status', 'unknown')
        logger.info(f"Null distribution computation status after join: {final_status}")
        
        if final_status == 'error':
            raise RuntimeError(f"Null distribution computation failed: {null_dist_result.get('error', 'Unknown error')}")
        elif final_status == 'unavailable':
            print(f"\n⚠ Statistical benchmarking unavailable: {null_dist_result.get('error', 'Unknown reason')}")
            print("   Continuing with enrichment analysis only (no statistical benchmarking)")
            null_distribution = None
        elif final_status == 'running':
            logger.warning("Null distribution computation still running after timeout, proceeding anyway...")
            null_distribution = None
        elif final_status == 'completed':
            # Get null distribution (convert to expected format)
            if null_dist_result.get('null_distribution'):
                null_distribution = null_dist_result['null_distribution']
                permutation_p_threshold = null_dist_result.get('permutation_p_threshold', 0.05)
                user_p_threshold = null_dist_result.get('user_p_threshold', PARAMS["p_threshold"])
                
                # Check if null_distribution is actually populated
                if null_distribution and len(null_distribution) > 0:
                    print(f"✓ Null distribution computed for {len(null_distribution)} gene list size(s)")
                    if user_p_threshold != permutation_p_threshold:
                        print(f"  Note: Permutation data filtered to p-value <= {user_p_threshold}")
                        print(f"        (Original permutation data generated with p-value threshold {permutation_p_threshold})")
                else:
                    logger.warning(f"Null distribution computation completed but result is empty. Status: {final_status}, keys: {list(null_dist_result.keys())}")
                    null_distribution = None
            else:
                logger.warning(f"Null distribution computation completed but no result found. Status: {final_status}, keys: {list(null_dist_result.keys())}")
                null_distribution = None
        else:
            logger.warning(f"Unexpected status: {final_status}, keys: {list(null_dist_result.keys())}")
            null_distribution = None
    elif null_dist_result['status'] == 'unavailable':
        print(f"\n⚠ Statistical benchmarking unavailable: {null_dist_result.get('error', 'Unknown reason')}")
        print("   Continuing with enrichment analysis only (no statistical benchmarking)")
        null_distribution = None
    else:
        # Status is not 'running' or 'unavailable' - might be 'completed' already
        final_status = null_dist_result.get('status', 'unknown')
        logger.info(f"Null distribution status (not running): {final_status}")
        if final_status == 'completed' and null_dist_result.get('null_distribution'):
            null_distribution = null_dist_result['null_distribution']
            if null_distribution and len(null_distribution) > 0:
                print(f"✓ Null distribution already computed for {len(null_distribution)} gene list size(s)")
            else:
                null_distribution = None
        else:
            null_distribution = None
    
    # 7. Combined network analysis (FULL network with all libraries)
    print("\n" + "=" * 80)
    print("COMBINED NETWORK ANALYSIS (All Libraries)")
    print("=" * 80)
    print()
    print("Note: Individual libraries are not analyzed separately because")
    print("      by design, terms in a single library are not connected.")
    print("      Only the combined network across all libraries has meaningful clusters.")
    print()
    
    # Combine ALL iGEA results (full network for user)
    combined_analyzer = NetworkConnectivityAnalyzer()
    for lib_name, iter_enrich in all_results.items():
        if iter_enrich is None:
            continue
        
        results = []
        for record in iter_enrich.results:
            genes_removed = record.get("Genes removed for next iteration", [])
            if isinstance(genes_removed, str):
                genes_removed = [g.strip() for g in genes_removed.split(",") if g.strip()]
            
            results.append({
                'Term': f"{lib_name}: {record.get('Term', '')}",
                'Iteration': record.get("Iteration", 1),
                'Library': lib_name,  # Add library name for diversity tracking
                'Genes removed for next iteration': genes_removed,
            })
        
        combined_analyzer.add_igea_results(results)
    
    # Get clusters from FULL network (sorted by size, largest first)
    clusters_full = combined_analyzer.get_clusters()
    
    print(f"Found {len(clusters_full)} clusters in FULL combined network (all {len(libraries)} libraries)")
    print()
    
    # 8. Create filtered network for statistics (only libraries with permutation data)
    if null_distribution is None:
        print("\n" + "=" * 80)
        print("STATISTICAL BENCHMARKING")
        print("=" * 80)
        print()
        print("⚠ Statistical benchmarking is not available.")
        print(f"   Reason: {null_dist_result.get('error', 'Unknown')}")
        print()
        print("   Continuing with network analysis only (no statistical comparison)")
        print()
        # Skip benchmarking sections
        filtered_benchmark = None
        clusters_for_stats = []
    else:
        print("=" * 80)
        print("STATISTICAL BENCHMARKING")
        print("=" * 80)
        print()
        print(f"⚠ IMPORTANT: Statistics computed using only {len(libraries_with_data)} libraries with permutation data:")
        print(f"   {', '.join(libraries_with_data)}")
        if libraries_without_data:
            print(f"\n   Libraries included in enrichment but EXCLUDED from statistics:")
            print(f"   {', '.join(libraries_without_data)}")
            print(f"   (Permutation data not available for these libraries)")
        print()
    
    # Create filtered analyzer with only libraries that have permutation data
    # Need to match user library names (e.g., "Reactome") to parquet library names (e.g., "C2: CP: Reactome Pathways")
    filtered_analyzer = NetworkConnectivityAnalyzer()
    
    # Create mapping from user library names to parquet library names
    user_to_parquet_mapping = {}
    for user_lib in all_results.keys():
        normalized_user = normalize_library_name(user_lib)
        for parquet_lib in libraries_with_data:
            normalized_parquet = normalize_library_name(parquet_lib)
            if (normalized_user == normalized_parquet or 
                normalized_user in normalized_parquet.lower() or 
                normalized_parquet in normalized_user.lower()):
                user_to_parquet_mapping[user_lib] = parquet_lib
                break
    
    for lib_name, iter_enrich in all_results.items():
        if iter_enrich is None:
            continue
        # Only include libraries that match libraries with permutation data
        if lib_name not in user_to_parquet_mapping:
            continue
        
        results = []
        for record in iter_enrich.results:
            genes_removed = record.get("Genes removed for next iteration", [])
            if isinstance(genes_removed, str):
                genes_removed = [g.strip() for g in genes_removed.split(",") if g.strip()]
            
            results.append({
                'Term': f"{lib_name}: {record.get('Term', '')}",
                'Iteration': record.get("Iteration", 1),
                'Library': lib_name,
                'Genes removed for next iteration': genes_removed,
            })
        
        filtered_analyzer.add_igea_results(results)
    
    # Get clusters from filtered network (for statistics)
    clusters_filtered = filtered_analyzer.get_clusters()
    
    print(f"Filtered network: {len(clusters_filtered)} clusters (using {len(libraries_with_data)} libraries)")
    print()
    
    # Convert filtered results to benchmark format
    filtered_results = []
    for lib_name, iter_enrich in all_results.items():
        if iter_enrich is None or lib_name not in user_to_parquet_mapping:
            continue
        for record in iter_enrich.results:
            genes_removed = record.get("Genes removed for next iteration", [])
            if isinstance(genes_removed, str):
                genes_removed = [g.strip() for g in genes_removed.split(",") if g.strip()]
            filtered_results.append({
                'Term': f"{lib_name}: {record.get('Term', '')}",
                'Iteration': record.get("Iteration", 1),
                'Library': lib_name,
                'Genes removed for next iteration': genes_removed,
            })
    
    # Get filtered network metrics for benchmarking
    if null_distribution is not None:
        filtered_metrics = filtered_analyzer.compute_metrics(include_library_diversity=True)
        filtered_benchmark = benchmark_real_results(
            filtered_results,
            gene_set.size,
            null_distribution,
            use_interpolation=True
        )
        
        # Print filtered network benchmark
        print("=" * 80)
        print("NETWORK BENCHMARK (Statistics Libraries Only)")
        print("=" * 80)
        print_benchmark_results(filtered_benchmark, f"FILTERED NETWORK ({len(libraries_with_data)} libraries)")
    else:
        filtered_metrics = filtered_analyzer.compute_metrics(include_library_diversity=True)
        filtered_benchmark = None
    
    # 9. Generate cluster-by-cluster statistical report
    if null_distribution is not None:
        # Use FULL network clusters but filter to only those containing statistics libraries
        print("\n" + "=" * 80)
        print("CLUSTER-BY-CLUSTER STATISTICAL REPORT")
        print("=" * 80)
        print()
        print(f"Note: Statistics computed using {len(libraries_with_data)} libraries: {', '.join(libraries_with_data)}")
        if libraries_without_data:
            print(f"      Full network includes {len(libraries)} libraries total (including {', '.join(libraries_without_data)})")
        print()
        
        # Import benchmark function
        from code.network_connectivity_benchmark import benchmark_cluster
        
        # Benchmark only the largest cluster (null distribution is built from largest clusters only)
        # Clusters are already sorted by size (largest first)
        # Filter to only the largest cluster that contains at least one library with permutation data
        clusters_for_stats = []
        
        # Find the largest cluster that has matching libraries
        for cluster in clusters_full:
            # Check if cluster contains any library with permutation data
            libraries_in_cluster = set()
            for term in cluster['terms']:
                library = combined_analyzer.term_to_library.get(term, "Unknown")
                libraries_in_cluster.add(library)
            
            # Check if any library in cluster matches libraries with permutation data
            # Use normalized matching like in benchmarking.py
            cluster_has_matching_lib = False
            for cluster_lib in libraries_in_cluster:
                normalized_cluster = normalize_library_name(cluster_lib)
                for parquet_lib in libraries_with_data:
                    normalized_parquet = normalize_library_name(parquet_lib)
                    if (normalized_cluster == normalized_parquet or 
                        normalized_cluster in normalized_parquet.lower() or 
                        normalized_parquet in normalized_cluster.lower()):
                        cluster_has_matching_lib = True
                        break
                if cluster_has_matching_lib:
                    break
            
            if cluster_has_matching_lib:
                clusters_for_stats.append(cluster)
                # Only benchmark the largest cluster (first one that matches)
                break
        
        if clusters_for_stats:
            print(f"Benchmarking largest cluster: {clusters_for_stats[0]['size']} nodes "
                  f"({clusters_for_stats[0]['n_genes']} genes, {clusters_for_stats[0]['n_terms']} terms)")
            print(f"Note: Only the largest cluster is benchmarked against null distribution "
                  f"(which is built from largest clusters in each permutation)")
        else:
            print(f"No clusters found that contain statistics libraries")
        print()
    else:
        clusters_for_stats = []
    
    # Build cluster statistics table with benchmark results (only if null distribution available)
    cluster_rows = []
    if null_distribution is not None and clusters_for_stats:
        for cluster in clusters_for_stats:
            # Get terms in this cluster (with library info and centrality)
            term_data = []
            libraries_in_cluster = set()
            for term in cluster['terms']:
                library = combined_analyzer.term_to_library.get(term, "Unknown")
                libraries_in_cluster.add(library)
                # Extract term name (remove library prefix if present)
                term_name = term.split(": ", 1)[-1] if ": " in term else term
                # Calculate term centrality (degree = number of genes connected to this term)
                centrality = len(combined_analyzer.term_to_genes.get(term, set()))
                term_data.append({
                    'term': f"{library}:{term_name}",
                    'centrality': centrality,
                    'original_term': term
                })
            
            # Rank terms by centrality (highest first)
            term_data.sort(key=lambda x: x['centrality'], reverse=True)
            term_list = [t['term'] for t in term_data]
            term_centralities = {t['term']: t['centrality'] for t in term_data}
            
            # Benchmark this cluster
            # Debug: Check null_distribution format
            if not null_distribution:
                logger.warning(f"null_distribution is empty or None when calling benchmark_cluster")
            else:
                logger.info(f"null_distribution keys: {list(null_distribution.keys())}")
                for size_key, stats in null_distribution.items():
                    if isinstance(stats, dict):
                        logger.info(f"  Size {size_key}: {len(stats)} metrics")
                        logger.info(f"    Metric keys: {list(stats.keys())[:10]}")
            
            cluster_benchmark = benchmark_cluster(
                cluster,
                null_distribution,
                gene_set.size,
                use_interpolation=True
            )
            
            # Debug: Check what was returned
            if not cluster_benchmark:
                logger.warning(f"benchmark_cluster returned empty result")
            else:
                logger.info(f"benchmark_cluster returned {len(cluster_benchmark)} metrics: {list(cluster_benchmark.keys())}")
            
            # Extract benchmark statistics (including null distribution stats)
            size_bench = cluster_benchmark.get('cluster_size', {})
            genes_bench = cluster_benchmark.get('cluster_genes', {})
            terms_bench = cluster_benchmark.get('cluster_terms', {})
            edges_bench = cluster_benchmark.get('cluster_edges', {})
            density_bench = cluster_benchmark.get('cluster_avg_edges_per_gene', cluster_benchmark.get('cluster_density', {}))
            lib_bench = cluster_benchmark.get('cluster_libraries', {})
            
            cluster_rows.append({
                'Cluster_Number': cluster['cluster_number'],
                'Cluster_Size': cluster['size'],
                'Cluster_Size_Z': size_bench.get('z_score', 0.0),
                'Cluster_Size_Pct': size_bench.get('percentile', 50.0),
                'Cluster_Size_Mean': size_bench.get('null_mean', 0.0),
                'Cluster_Size_Std': size_bench.get('null_std', 0.0),
                'N_Genes': cluster['n_genes'],
                'N_Genes_Z': genes_bench.get('z_score', 0.0),
                'N_Genes_Pct': genes_bench.get('percentile', 50.0),
                'N_Genes_Mean': genes_bench.get('null_mean', 0.0),
                'N_Genes_Std': genes_bench.get('null_std', 0.0),
                'N_Terms': cluster['n_terms'],
                'N_Terms_Z': terms_bench.get('z_score', 0.0),
                'N_Terms_Pct': terms_bench.get('percentile', 50.0),
                'N_Terms_Mean': terms_bench.get('null_mean', 0.0),
                'N_Terms_Std': terms_bench.get('null_std', 0.0),
                'N_Edges': cluster['n_edges'],
                'N_Edges_Z': edges_bench.get('z_score', 0.0),
                'N_Edges_Pct': edges_bench.get('percentile', 50.0),
                'N_Edges_Mean': edges_bench.get('null_mean', 0.0),
                'N_Edges_Std': edges_bench.get('null_std', 0.0),
                'Average_Edges_per_Gene': cluster.get('avg_edges_per_gene', cluster.get('density', 0)),
                'Average_Edges_per_Gene_Z': density_bench.get('z_score', 0.0),
                'Average_Edges_per_Gene_Pct': density_bench.get('percentile', 50.0),
                'Average_Edges_per_Gene_Mean': density_bench.get('null_mean', 0.0),
                'Average_Edges_per_Gene_Std': density_bench.get('null_std', 0.0),
                # Keep old column names for backward compatibility
                'Cluster_Density': cluster.get('avg_edges_per_gene', cluster.get('density', 0)),
                'Cluster_Density_Z': density_bench.get('z_score', 0.0),
                'Cluster_Density_Pct': density_bench.get('percentile', 50.0),
                'Cluster_Density_Mean': density_bench.get('null_mean', 0.0),
                'Cluster_Density_Std': density_bench.get('null_std', 0.0),
                'N_Libraries': len(libraries_in_cluster),
                'N_Libraries_Z': lib_bench.get('z_score', 0.0),
                'N_Libraries_Pct': lib_bench.get('percentile', 50.0),
                'N_Libraries_Mean': lib_bench.get('null_mean', 0.0),
                'N_Libraries_Std': lib_bench.get('null_std', 0.0),
                'Libraries': ', '.join(sorted(libraries_in_cluster)),
                'Terms': '; '.join(term_list),
                'Term_Centralities': '; '.join([f"{t['term']}:{t['centrality']}" for t in term_data]),
                'N_Terms_Total': len(term_list),
            })
    
    # Create DataFrame (only if we have cluster rows)
    if cluster_rows:
        cluster_df = pd.DataFrame(cluster_rows)
    else:
        cluster_df = pd.DataFrame()  # Empty DataFrame
    
    # Print cluster table with statistics (only if we have data)
    if len(cluster_df) > 0:
        print("Largest Cluster Statistical Benchmarks:")
        print("(Only the largest cluster is benchmarked - null distribution built from largest clusters only)")
        print("=" * 120)
        print()
        
        # Format for display
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', 50)
        
        # Print summary table (key metrics)
        summary_cols = ['Cluster_Number', 'Cluster_Size', 'Cluster_Size_Z', 'Cluster_Size_Pct',
                        'N_Genes', 'N_Genes_Z', 'N_Terms', 'N_Terms_Z',
                        'Cluster_Density', 'Cluster_Density_Z', 'N_Libraries', 'N_Libraries_Z']
        summary_df = cluster_df[summary_cols].copy()
        
        # Format columns for better readability
        for col in summary_df.columns:
            if col.endswith('_Z'):
                summary_df[col] = summary_df[col].apply(lambda x: f"{x:>6.2f}")
            elif col.endswith('_Pct'):
                summary_df[col] = summary_df[col].apply(lambda x: f"{x:>5.1f}%")
            elif 'Density' in col:
                summary_df[col] = summary_df[col].apply(lambda x: f"{x:>6.4f}")
            else:
                summary_df[col] = summary_df[col].apply(lambda x: f"{x:>6.0f}" if pd.notna(x) and isinstance(x, (int, float)) else str(x))
        
        print(summary_df.to_string(index=False))
        print()
        
        # Print full table
        print("=" * 120)
        print("Full Cluster Statistics:")
        print("-" * 120)
        print(cluster_df.to_string(index=False))
        print()
    else:
        print("No cluster statistics available (statistical benchmarking not performed)")
        print()
    
    # 9. Save results
    print("=" * 80)
    print("Saving results...")
    print("=" * 80)
    
    output_dir = PROJECT_ROOT / "results" / "hiv_connectivity_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize cluster_file variable
    cluster_file = None
    
    # Save cluster-by-cluster results (with statistics) if available
    if len(cluster_df) > 0:
        cluster_file = output_dir / "hiv_clusters_by_size.tsv"
        cluster_df.to_csv(cluster_file, sep="\t", index=False)
        print(f"✓ Saved cluster-by-cluster statistical report to: {cluster_file}")
        
        # Also save a summary version (key metrics only)
        summary_file = output_dir / "hiv_clusters_summary.tsv"
        summary_cols = ['Cluster_Number', 'Cluster_Size', 'Cluster_Size_Z', 'Cluster_Size_Pct',
                        'N_Genes', 'N_Genes_Z', 'N_Terms', 'N_Terms_Z',
                        'Cluster_Density', 'Cluster_Density_Z', 'N_Libraries', 'N_Libraries_Z',
                        'Libraries']
        summary_df = cluster_df[summary_cols].copy()
        summary_df.to_csv(summary_file, sep="\t", index=False)
        print(f"✓ Saved cluster summary (key metrics) to: {summary_file}")
        
        # Generate formatted statistical report
        try:
            from generate_cluster_statistical_report import generate_cluster_report
            report_file = output_dir / "hiv_clusters_statistical_report.txt"
            generate_cluster_report(str(cluster_file), str(report_file), gene_list_name="HIV Gene List")
            print(f"✓ Generated formatted statistical report: {report_file}")
            
            # Add library information note to report
            with open(report_file, 'r') as f:
                report_content = f.read()
            
            # Add note about libraries used for statistics
            library_note = f"""
====================================================================================================
IMPORTANT: Library Information for Statistical Analysis
====================================================================================================

Statistics were computed using {len(libraries_with_data)} libraries with permutation data:
  {', '.join(libraries_with_data)}

"""
            if libraries_without_data:
                library_note += f"""Libraries included in enrichment but EXCLUDED from statistics:
  {', '.join(libraries_without_data)}
  (Permutation data not available for these libraries)

"""
            library_note += """The full network visualization includes all selected libraries.
However, statistical benchmarks are only computed for libraries with available
permutation data to ensure accurate comparison against the null distribution.

====================================================================================================

"""
            
            # Insert note after the header
            lines = report_content.split('\n')
            insert_pos = 0
            for i, line in enumerate(lines):
                if 'Gene List:' in line or 'Total Clusters:' in line:
                    insert_pos = i + 3
                    break
            
            # Insert library note
            for i, note_line in enumerate(library_note.strip().split('\n')):
                lines.insert(insert_pos + i, note_line)
            
            with open(report_file, 'w') as f:
                f.write('\n'.join(lines))
            
            print(f"✓ Added library information note to report")
        except Exception as e:
            logger.warning(f"Could not generate formatted report: {e}")
            import traceback
            traceback.print_exc()
    
    # Save overall benchmark results
    benchmark_file = output_dir / "hiv_connectivity_benchmark.json"
    benchmark_data = {
        "gene_list_size": gene_set.size,
        "libraries_selected": user_selected_libraries,
        "libraries_with_permutation_data": libraries_with_data,
        "libraries_without_permutation_data": libraries_without_data,
        "note": "Full network includes all selected libraries. Statistics computed using only libraries with permutation data.",
        "full_network": {
            "n_clusters": len(clusters_full),
            "metrics": combined_analyzer.compute_metrics(include_library_diversity=True),
        },
        "filtered_network": {
            "n_clusters": len(clusters_filtered),
            "metrics": filtered_metrics,
            "benchmark": filtered_benchmark if filtered_benchmark is not None else "Not available (p-value threshold > 0.05)",
        },
        "clusters": cluster_rows,
    }
    with open(benchmark_file, 'w') as f:
        json.dump(benchmark_data, f, indent=2, default=str)
    print(f"✓ Saved benchmark results to: {benchmark_file}")
    
    print("\n" + "=" * 80)
    print("Benchmarking complete!")
    print("=" * 80)
    print()
    print(f"Summary:")
    print(f"  Full network: {len(clusters_full)} clusters (all {len(libraries)} libraries)")
    if null_distribution is not None:
        print(f"  Statistics network: {len(clusters_filtered)} clusters ({len(libraries_with_data)} libraries with permutation data)")
        if clusters_for_stats:
            print(f"  Largest cluster: {clusters_for_stats[0]['size']} nodes ({clusters_for_stats[0]['n_genes']} genes, {clusters_for_stats[0]['n_terms']} terms)")
    else:
        print(f"  Statistical benchmarking: Not available (p-value threshold > 0.05)")
        print(f"  Filtered network: {len(clusters_filtered)} clusters ({len(libraries_with_data)} libraries)")
        if clusters_full:
            print(f"  Largest cluster: {clusters_full[0]['size']} nodes ({clusters_full[0]['n_genes']} genes, {clusters_full[0]['n_terms']} terms)")
    
    if cluster_file and cluster_file.exists():
        print(f"         See {cluster_file} for detailed cluster-by-cluster statistics")


if __name__ == "__main__":
    main()
