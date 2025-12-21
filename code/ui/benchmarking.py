"""
Benchmarking functionality for Streamlit app.
Computes network connectivity benchmarks and displays statistical tables.
"""

import logging
import threading
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd

from network_connectivity_benchmark import (
    NetworkConnectivityAnalyzer,
    benchmark_cluster
)
from parallel_null_distribution import (
    find_intersection_libraries,
    compute_null_distribution_parallel,
    normalize_library_name
)

logger = logging.getLogger(__name__)


def format_status(z_score: float, percentile: float, is_available: bool = True) -> str:
    """
    Format status indicator based on z-score and percentile.
    
    Args:
        z_score: Z-score value
        percentile: Percentile value
        is_available: Whether the metric is available in null distribution
        
    Returns:
        Status string
    """
    if not is_available or (z_score == 0.0 and percentile == 50.0):
        return "Not available"
    
    if z_score > 2.0 and percentile > 95.0:
        return "✓✓ SIGNIFICANTLY BETTER"
    elif z_score > 1.0 and percentile > 84.0:
        return "✓ BETTER"
    elif z_score > -1.0 and percentile > 16.0:
        return "~ SIMILAR"
    elif z_score < -1.0 and percentile < 16.0:
        return "✗ WORSE"
    else:
        return "~ SIMILAR"


def compute_benchmark_for_streamlit(
    iter_enrich: Dict,
    gene_list_size: int,
    p_threshold: float,
    parquet_dir: Path,
    merged_permutation_file: Optional[Path] = None
) -> Tuple[Optional[Dict], Optional[List[Dict]], List[str], List[str]]:
    """
    Compute benchmarking for Streamlit app.
    
    Args:
        iter_enrich: Dictionary of IterativeEnrichment objects by library name
        gene_list_size: Size of the input gene list
        p_threshold: P-value threshold used for enrichment
        parquet_dir: Directory containing Parquet files with permutation statistics
        merged_permutation_file: Optional path to merged permutation TSV file
        
    Returns:
        Tuple of (null_distribution, cluster_benchmarks, libraries_with_data, libraries_without_data, actual_size_used)
        Returns None for null_distribution and cluster_benchmarks if benchmarking unavailable
        actual_size_used: The gene list size from permutation data that was actually used
    """
    # Get all library names from iter_enrich
    user_selected_libraries = list(iter_enrich.keys())
    
    # Find which libraries have permutation data
    libraries_with_data, libraries_without_data = find_intersection_libraries(
        user_selected_libraries,
        parquet_dir
    )
    
    if not libraries_with_data:
        logger.warning("No libraries with permutation data available for benchmarking")
        return None, None, [], user_selected_libraries, gene_list_size
    
    # Check if p-value threshold is valid for benchmarking
    if p_threshold > 0.05:
        logger.warning(f"P-value threshold {p_threshold} > 0.05, benchmarking may be unavailable")
        # Still try to compute, but warn user
    
    # Compute null distribution in parallel
    null_dist_result = {
        'null_distribution': None,
        'status': 'running',
        'libraries_used': libraries_with_data,
        'error': None
    }
    null_dist_lock = threading.Lock()
    
    # Try to find merged permutation file
    if merged_permutation_file is None:
        # Try multiple possible locations
        possible_paths = [
            parquet_dir.parent / "permutation_results" / "merged_permutation_results.tsv",
            parquet_dir.parent.parent / "archive" / "permutation_analysis" / "results" / "permutation_results" / "merged_permutation_results.tsv",
            parquet_dir.parent.parent / "merged_permutation_results.tsv",
            parquet_dir.parent.parent.parent / "merged_permutation_results.tsv",  # In case parquet_dir is nested deeper
        ]
        for path in possible_paths:
            if path.exists():
                merged_permutation_file = path
                logger.info(f"Found merged permutation file: {merged_permutation_file}")
                break
        else:
            logger.warning("Could not find merged permutation file, proceeding without p-value filtering")
    
    null_dist_thread = threading.Thread(
        target=compute_null_distribution_parallel,
        args=(
            parquet_dir,
            gene_list_size,
            libraries_with_data,
            null_dist_result,
            null_dist_lock,
            p_threshold,
            merged_permutation_file
        ),
        daemon=True
    )
    null_dist_thread.start()
    
    # Wait for null distribution computation
    null_dist_thread.join(timeout=120)  # Wait up to 2 minutes
    
    # Check status
    final_status = null_dist_result.get('status', 'unknown')
    if final_status == 'error':
        logger.error(f"Null distribution computation failed: {null_dist_result.get('error', 'Unknown error')}")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size
    elif final_status == 'unavailable':
        logger.warning(f"Statistical benchmarking unavailable: {null_dist_result.get('error', 'Unknown reason')}")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size
    elif final_status == 'running':
        logger.warning("Null distribution computation still running after timeout")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size
    elif final_status == 'completed':
        null_distribution = null_dist_result.get('null_distribution')
        if not null_distribution or len(null_distribution) == 0:
            logger.warning("Null distribution computation completed but result is empty")
            return None, None, libraries_with_data, libraries_without_data, gene_list_size
    else:
        logger.warning(f"Unexpected status: {final_status}")
        return None, None, libraries_with_data, libraries_without_data, gene_list_size
    
    # Build combined network from ONLY libraries with permutation data
    # This ensures fair comparison: real network uses same libraries as null distribution
    # We filter iGEA results to match the libraries used in the null distribution
    combined_analyzer = NetworkConnectivityAnalyzer()
    
    # Create a mapping from normalized library names to original names for matching
    # Map user library names (from iter_enrich) to permutation library names (from libraries_with_data)
    lib_name_mapping = {}
    for lib_name in iter_enrich.keys():
        normalized = normalize_library_name(lib_name)
        lib_name_mapping[lib_name] = normalized
    
    # Find which libraries from iter_enrich match libraries_with_data
    libraries_to_use = []
    for lib_name in iter_enrich.keys():
        normalized_user = lib_name_mapping[lib_name]
        # Check if this library matches any library with permutation data
        for lib_with_data in libraries_with_data:
            normalized_data = normalize_library_name(lib_with_data)
            # Match if normalized names are the same or one contains the other
            if (normalized_user == normalized_data or 
                normalized_user in normalized_data.lower() or 
                normalized_data in normalized_user.lower()):
                if lib_name not in libraries_to_use:
                    libraries_to_use.append(lib_name)
                    break
    
    if not libraries_to_use:
        logger.warning(f"Could not match any libraries from iter_enrich to libraries with permutation data. "
                      f"User libraries: {list(iter_enrich.keys())}, "
                      f"Permutation libraries: {libraries_with_data}")
        # Fallback: try to use all libraries_with_data (but this might not work if names don't match)
        libraries_to_use = list(iter_enrich.keys())
    
    logger.info(f"Filtering iGEA results: Using {len(libraries_to_use)} libraries for benchmarking "
               f"(matching {len(libraries_with_data)} libraries with permutation data): {libraries_to_use}")
    
    for lib_name in libraries_to_use:
        iter_enrich_obj = iter_enrich.get(lib_name)
        if iter_enrich_obj is None:
            continue
        
        results = []
        for record in iter_enrich_obj.results:
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
    
    # Get clusters
    clusters = combined_analyzer.get_clusters()
    
    if not clusters:
        logger.info("No clusters found in network")
        # Get actual size used
        actual_size_used = gene_list_size
        if null_distribution:
            available_sizes = sorted([int(k) for k in null_distribution.keys()])
            if available_sizes:
                actual_size_used = min(available_sizes, key=lambda x: abs(x - gene_list_size))
        return null_distribution, [], libraries_with_data, libraries_without_data, actual_size_used
    
    # Benchmark each cluster
    cluster_benchmarks = []
    for cluster in clusters:
        benchmark = benchmark_cluster(
            cluster,
            null_distribution,
            gene_list_size,
            use_interpolation=True
        )
        
        if benchmark:
            cluster_benchmarks.append({
                'cluster': cluster,
                'benchmark': benchmark
            })
    
    # Get the actual size used from null distribution
    actual_size_used = gene_list_size
    if null_distribution:
        available_sizes = sorted([int(k) for k in null_distribution.keys()])
        if available_sizes:
            actual_size_used = min(available_sizes, key=lambda x: abs(x - gene_list_size))
    
    return null_distribution, cluster_benchmarks, libraries_with_data, libraries_without_data, actual_size_used


def extract_benchmark_table_data(cluster_benchmarks: List[Dict]) -> List[Dict]:
    """
    Extract table data for "Statistical Benchmarks vs Random Gene Lists" table.
    
    Args:
        cluster_benchmarks: List of cluster benchmark dictionaries
        
    Returns:
        List of dictionaries with table row data (one per cluster, showing largest cluster only)
    """
    if not cluster_benchmarks:
        return []
    
    # For now, show only the largest cluster (first one, as they're sorted by size)
    # In the future, we could show all clusters or allow user to select
    largest_cluster_data = cluster_benchmarks[0]
    cluster = largest_cluster_data['cluster']
    benchmark = largest_cluster_data['benchmark']
    
    table_rows = []
    
    # Define metrics to display
    metrics = [
        ('Cluster Size', 'cluster_size', int),
        ('Number of Genes', 'cluster_genes', int),
        ('Number of Terms', 'cluster_terms', int),
        ('Number of Edges', 'cluster_edges', int),
        ('Average Edges per Gene', 'cluster_avg_edges_per_gene', float),
        ('Number of Libraries', 'cluster_libraries', int),
    ]
    
    for metric_name, metric_key, value_type in metrics:
        if metric_key not in benchmark:
            continue
        
        metric_data = benchmark[metric_key]
        real_value = metric_data['real_value']
        null_mean = metric_data.get('null_mean', 0.0)
        null_std = metric_data.get('null_std', 0.0)
        z_score = metric_data.get('z_score', 0.0)
        percentile = metric_data.get('percentile', 50.0)
        
        # Check if metric is available
        is_available = not (z_score == 0.0 and percentile == 50.0 and null_mean == 0.0)
        status = format_status(z_score, percentile, is_available)
        
        # Format value based on type
        if value_type == int:
            value_str = f"{int(real_value)}"
            mean_str = f"{null_mean:.2f}" if is_available and null_mean > 0 else "N/A"
            std_str = f"{null_std:.2f}" if is_available and null_std > 0 else "N/A"
        else:  # float
            value_str = f"{real_value:.4f}"
            mean_str = f"{null_mean:.4f}" if is_available and null_mean > 0 else "N/A"
            std_str = f"{null_std:.4f}" if is_available and null_std > 0 else "N/A"
        
        table_rows.append({
            'Metric': metric_name,
            'Value': value_str,
            'Null Mean': mean_str,
            'Null Std': std_str,
            'Z-score': f"{z_score:.2f}",
            'Percentile': f"{percentile:.1f}%",
            'Status': status
        })
    
    return table_rows


def generate_statistical_report_text(
    cluster_benchmarks: List[Dict],
    gene_list_name: str,
    libraries_with_data: List[str],
    libraries_without_data: List[str]
) -> str:
    """
    Generate full statistical report text (same format as generate_cluster_statistical_report.py).
    
    Args:
        cluster_benchmarks: List of cluster benchmark dictionaries
        gene_list_name: Name of the gene list
        libraries_with_data: Libraries used for statistics
        libraries_without_data: Libraries excluded from statistics
        
    Returns:
        Full report text as string
    """
    lines = []
    
    # Header
    lines.append("=" * 100)
    lines.append("CLUSTER-BY-CLUSTER STATISTICAL REPORT")
    lines.append("=" * 100)
    lines.append("")
    lines.append(f"Gene List: {gene_list_name}")
    lines.append(f"Total Clusters: {len(cluster_benchmarks)}")
    lines.append("")
    
    # Library information
    lines.append("=" * 100)
    lines.append("IMPORTANT: Library Information for Statistical Analysis")
    lines.append("=" * 100)
    lines.append("")
    lines.append(f"Statistics were computed using {len(libraries_with_data)} libraries with permutation data:")
    lines.append(f"  {', '.join(libraries_with_data)}")
    lines.append("")
    if libraries_without_data:
        lines.append(f"Libraries included in enrichment but EXCLUDED from statistics:")
        lines.append(f"  {', '.join(libraries_without_data)}")
        lines.append(f"  (Permutation data not available for these libraries)")
        lines.append("")
    lines.append("The full network visualization includes all selected libraries.")
    lines.append("However, statistical benchmarks are only computed for libraries with available")
    lines.append("permutation data to ensure accurate comparison against the null distribution.")
    lines.append("")
    lines.append("=" * 100)
    lines.append("")
    lines.append("This report shows each cluster with benchmark statistics comparing")
    lines.append("against random gene lists of similar size.")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("  • Z-score > 2.0 and Percentile > 95%: Significantly better than random (top 5%)")
    lines.append("  • Z-score > 1.0 and Percentile > 84%: Better than random (top 16%)")
    lines.append("  • Z-score ~ 0.0 and Percentile ~ 50%: Similar to random")
    lines.append("  • Z-score < -1.0 and Percentile < 16%: Worse than random")
    lines.append("")
    lines.append("=" * 100)
    lines.append("")
    
    # Process each cluster
    for idx, cluster_data in enumerate(cluster_benchmarks, 1):
        cluster = cluster_data['cluster']
        benchmark = cluster_data['benchmark']
        
        lines.append("-" * 100)
        lines.append(f"CLUSTER {idx} (Ranked by Size)")
        lines.append("-" * 100)
        lines.append("")
        
        # Basic metrics
        lines.append("Basic Cluster Metrics:")
        lines.append(f"  Cluster Size:        {cluster['size']} nodes (genes + terms)")
        lines.append(f"  Number of Genes:     {cluster['n_genes']}")
        lines.append(f"  Number of Terms:     {cluster['n_terms']}")
        lines.append(f"  Number of Edges:     {cluster['n_edges']}")
        lines.append(f"  Cluster Density:     {cluster['density']:.4f}")
        lines.append(f"  Number of Libraries: {cluster.get('n_libraries', 0)}")
        if 'libraries' in cluster:
            lines.append(f"  Libraries:           {', '.join(cluster['libraries'])}")
        lines.append("")
        
        # Statistical benchmarks
        lines.append("Statistical Benchmarks vs Random Gene Lists:")
        lines.append("  Metric                    Value      Null Mean   Null Std    Z-score   Percentile   Status")
        lines.append("  " + "-" * 110)
        
        # Define metrics
        metrics = [
            ('Cluster Size', 'cluster_size', int),
            ('Number of Genes', 'cluster_genes', int),
            ('Number of Terms', 'cluster_terms', int),
            ('Number of Edges', 'cluster_edges', int),
            ('Average Edges per Gene', 'cluster_avg_edges_per_gene', float),
            ('Number of Libraries', 'cluster_libraries', int),
        ]
        
        for metric_name, metric_key, value_type in metrics:
            if metric_key not in benchmark:
                continue
            
            metric_data = benchmark[metric_key]
            real_value = metric_data['real_value']
            null_mean = metric_data.get('null_mean', 0.0)
            null_std = metric_data.get('null_std', 0.0)
            z_score = metric_data.get('z_score', 0.0)
            percentile = metric_data.get('percentile', 50.0)
            
            is_available = not (z_score == 0.0 and percentile == 50.0 and null_mean == 0.0)
            status = format_status(z_score, percentile, is_available)
            
            if value_type == int:
                if is_available and null_mean > 0:
                    lines.append(f"  {metric_name:<25} {int(real_value):>8}  {null_mean:>10.2f}  {null_std:>10.2f}  {z_score:>8.2f}  {percentile:>8.1f}%  {status}")
                else:
                    lines.append(f"  {metric_name:<25} {int(real_value):>8}  {'N/A':>10}  {'N/A':>10}  {z_score:>8.2f}  {percentile:>8.1f}%  {status}")
            else:  # float
                if is_available and null_mean > 0:
                    lines.append(f"  {metric_name:<25} {real_value:>8.4f}  {null_mean:>10.4f}  {null_std:>10.4f}  {z_score:>8.2f}  {percentile:>8.1f}%  {status}")
                else:
                    lines.append(f"  {metric_name:<25} {real_value:>8.4f}  {'N/A':>10}  {'N/A':>10}  {z_score:>8.2f}  {percentile:>8.1f}%  {status}")
        
        lines.append("")
    
    # Footer
    lines.append("=" * 100)
    lines.append("END OF REPORT")
    lines.append("=" * 100)
    
    return '\n'.join(lines)

