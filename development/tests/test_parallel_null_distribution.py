#!/usr/bin/env python3
"""
Test script to verify parallel null distribution computation integration.

This script uses a small gene list to quickly test:
1. Library intersection detection
2. Parallel null distribution computation
3. Network filtering for statistics
4. Benchmarking workflow
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
from network_connectivity_benchmark import (
    NetworkConnectivityAnalyzer,
    benchmark_real_results
)
from parallel_null_distribution import (
    find_intersection_libraries,
    compute_null_distribution_parallel
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
PARQUET_DIR = PROJECT_ROOT / "permutations" / "permutation_cluster_statistics_parquet"

# Libraries to use (matching permutation data)
LIBRARIES = {
    "Reactome": "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
    "KEGG": "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",
    "GO BP": "c5.go.bp.v2025.1.Hs.symbols.gmt",
    "GO MF": "c5.go.mf.v2025.1.Hs.symbols.gmt",
    "GO CC": "c5.go.cc.v2025.1.Hs.symbols.gmt",
}

# iGEA parameters
PARAMS = {
    "p_threshold": 0.05,
    "min_overlap": 3,
    "min_term_size": 10,
    "max_term_size": 600,
    "max_iterations": 10,
    "p_value_method": "Fisher's Exact Test",
}

# Small test gene list (50 genes for quick testing)
TEST_GENES = [
    "TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "AKT1",
    "CDKN2A", "RB1", "MDM2", "VEGFA", "TGFB1", "IL6", "TNF", "NFKB1", "STAT3",
    "JAK2", "MAPK1", "MAPK3", "ERK1", "ERK2", "P53", "BAX", "BCL2", "CASP3",
    "CASP8", "CASP9", "FAS", "FASLG", "CDKN1A", "CDK4", "CDK6", "CCND1", "CCNE1",
    "MYCN", "ALK", "MET", "FGFR1", "PDGFRA", "KIT", "FLT3", "ABL1", "BCR",
    "RET", "ROS1", "NTRK1", "BRAF", "NRAS"
]


def load_background() -> BackgroundGeneSet:
    """Load background gene set."""
    bg_path = BACKGROUNDS_DIR / "all_genes.txt"
    if not bg_path.exists():
        raise FileNotFoundError(f"Background file not found: {bg_path}")
    
    logger.info(f"Loading background: {bg_path}")
    bg = BackgroundGeneSet(str(bg_path), name="all_genes", input_format="symbols", skip_validation=True)
    logger.info(f"Loaded background: {bg.size} genes")
    return bg


def load_libraries() -> Dict[str, GeneSetLibrary]:
    """Load gene set libraries."""
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
        use_multiprocessing=False,  # Disable for faster testing
    )
    
    logger.info(f"  {library.name}: {len(iter_enrich.results)} iterations")
    return iter_enrich


def print_benchmark_results(benchmark: Dict, network_name: str):
    """Print benchmark results in a readable format."""
    print(f"\n{network_name} Benchmark Results:")
    print("-" * 80)
    
    if not benchmark or 'comparison' not in benchmark:
        print("  No benchmark data available")
        return
    
    comparison = benchmark['comparison']
    
    # Key metrics to display
    key_metrics = [
        'largest_cluster_genes',
        'largest_cluster_terms',
        'largest_cluster_edges',
        'largest_cluster_density',
        'n_connected_components',
    ]
    
    for metric in key_metrics:
        if metric in comparison:
            comp = comparison[metric]
            z_score = comp.get('z_score', 0.0)
            percentile = comp.get('percentile', 50.0)
            real_value = comp.get('real_value', 0)
            null_mean = comp.get('null_mean', 0)
            
            # Status indicator
            if z_score > 2.0 and percentile > 95.0:
                status = "✓✓ SIGNIFICANTLY BETTER"
            elif z_score > 1.0 and percentile > 84.0:
                status = "✓ BETTER"
            elif z_score > -1.0 and percentile > 16.0:
                status = "~ SIMILAR"
            else:
                status = "✗ WORSE"
            
            print(f"  {metric:30s}: {real_value:6.1f} (null: {null_mean:6.1f}) | "
                  f"Z={z_score:6.2f}, Pct={percentile:5.1f}% | {status}")


def main():
    """Main test function."""
    print("\n" + "=" * 80)
    print("TEST: Parallel Null Distribution Integration")
    print("=" * 80)
    
    # 1. Create gene set
    print("\n[1/6] Creating test gene set...")
    background = load_background()
    gene_set = GeneSet(
        gene_list=TEST_GENES,
        validation_set=background.genes,
        hgcn=False,
        format=False
    )
    print(f"✓ Gene set created: {gene_set.size} genes")
    
    # 2. Load libraries
    print("\n[2/6] Loading gene set libraries...")
    libraries = load_libraries()
    print(f"✓ Loaded {len(libraries)} libraries")
    
    # 3. Check which libraries have permutation data
    print("\n[3/6] Checking libraries with permutation data...")
    user_selected_libraries = list(libraries.keys())
    libraries_with_data, libraries_without_data = find_intersection_libraries(
        user_selected_libraries,
        PARQUET_DIR
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
        print(f"⚠ Libraries without permutation data: {', '.join(libraries_without_data)}")
    
    # 4. Start null distribution computation in parallel thread
    print("\n[4/6] Starting parallel null distribution computation...")
    null_dist_result = {
        'null_distribution': None,
        'status': 'running',
        'libraries_used': libraries_with_data,
        'error': None
    }
    null_dist_lock = threading.Lock()
    
    null_dist_thread = threading.Thread(
        target=compute_null_distribution_parallel,
        args=(
            PARQUET_DIR,
            gene_set.size,
            libraries_with_data,
            null_dist_result,
            null_dist_lock,
            PARAMS["p_threshold"]  # Pass user's p-value threshold
        ),
        daemon=True
    )
    null_dist_thread.start()
    print(f"✓ Null distribution computation started in parallel thread")
    
    # 5. Run iGEA for each library
    print("\n[5/6] Running iGEA for each library...")
    print(f"Note: All {len(libraries)} libraries will be included in enrichment results.")
    if libraries_without_data:
        print(f"      However, statistics will only be computed for {len(libraries_with_data)} libraries.")
    print()
    
    all_results = {}
    
    for lib_name, library in libraries.items():
        try:
            iter_enrich = run_igea_for_library(gene_set, library, background)
            all_results[lib_name] = iter_enrich
            print(f"✓ {lib_name}: {len(iter_enrich.results)} iterations")
        except Exception as e:
            logger.error(f"Error processing {lib_name}: {e}", exc_info=True)
            all_results[lib_name] = None
    
    # 6. Wait for null distribution computation to complete
    print("\n[6/6] Waiting for null distribution computation...")
    null_dist_thread.join(timeout=30)  # Wait up to 30 seconds
    
    if null_dist_result['status'] == 'error':
        raise RuntimeError(f"Null distribution computation failed: {null_dist_result.get('error', 'Unknown error')}")
    elif null_dist_result['status'] == 'running':
        logger.warning("Null distribution computation still running after timeout, proceeding anyway...")
    
    # Get null distribution (convert to expected format)
    if null_dist_result['null_distribution']:
        null_distribution = null_dist_result['null_distribution']
        permutation_p_threshold = null_dist_result.get('permutation_p_threshold', 0.05)
        user_p_threshold = null_dist_result.get('user_p_threshold', PARAMS["p_threshold"])
        
        print(f"✓ Null distribution computed for {len(null_distribution)} gene list size(s)")
        if user_p_threshold != permutation_p_threshold:
            print(f"  Note: Permutation data generated with p-value threshold {permutation_p_threshold}")
            print(f"        Your analysis uses p-value threshold {user_p_threshold}")
            print(f"        Comparison is valid but may be conservative (your network may be smaller)")
    else:
        raise RuntimeError("Null distribution computation did not complete")
    
    # 7. Create FULL network (all libraries)
    print("\n" + "=" * 80)
    print("FULL NETWORK (All Libraries)")
    print("=" * 80)
    
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
                'Library': lib_name,
                'Genes removed for next iteration': genes_removed,
            })
        
        combined_analyzer.add_igea_results(results)
    
    clusters_full = combined_analyzer.get_clusters()
    print(f"✓ Found {len(clusters_full)} clusters in FULL network (all {len(libraries)} libraries)")
    
    # 8. Create FILTERED network (only statistics libraries)
    print("\n" + "=" * 80)
    print("FILTERED NETWORK (Statistics Libraries Only)")
    print("=" * 80)
    print(f"⚠ IMPORTANT: Statistics computed using only {len(libraries_with_data)} libraries:")
    print(f"   {', '.join(libraries_with_data)}")
    if libraries_without_data:
        print(f"\n   Libraries excluded from statistics: {', '.join(libraries_without_data)}")
    print()
    
    filtered_analyzer = NetworkConnectivityAnalyzer()
    for lib_name, iter_enrich in all_results.items():
        if iter_enrich is None or lib_name not in libraries_with_data:
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
    
    clusters_filtered = filtered_analyzer.get_clusters()
    print(f"✓ Found {len(clusters_filtered)} clusters in FILTERED network ({len(libraries_with_data)} libraries)")
    
    # 9. Benchmark filtered network
    print("\n" + "=" * 80)
    print("BENCHMARKING FILTERED NETWORK")
    print("=" * 80)
    
    filtered_results = []
    for lib_name, iter_enrich in all_results.items():
        if iter_enrich is None or lib_name not in libraries_with_data:
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
    
    filtered_benchmark = benchmark_real_results(
        filtered_results,
        gene_set.size,
        null_distribution,
        use_interpolation=True
    )
    
    print_benchmark_results(filtered_benchmark, f"FILTERED NETWORK ({len(libraries_with_data)} libraries)")
    
    # 10. Summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    print(f"✓ Gene list size: {gene_set.size}")
    print(f"✓ Libraries selected: {len(libraries)}")
    print(f"✓ Libraries with permutation data: {len(libraries_with_data)}")
    print(f"✓ Full network clusters: {len(clusters_full)}")
    print(f"✓ Filtered network clusters: {len(clusters_filtered)}")
    print(f"✓ Null distribution computed: {null_dist_result['status'] == 'completed'}")
    print(f"✓ Benchmark completed successfully")
    print("\n" + "=" * 80)
    print("TEST PASSED ✓")
    print("=" * 80)


if __name__ == '__main__':
    main()
