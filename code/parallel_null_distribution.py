#!/usr/bin/env python3
"""
Parallel null distribution computation from Parquet cluster statistics.

This module provides functions to:
1. Check which libraries have permutation data available
2. Compute null distribution from Parquet files in parallel with iGEA
3. Filter by user-selected libraries
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict, Set, Optional
from concurrent.futures import ThreadPoolExecutor
import threading

logger = logging.getLogger(__name__)


def get_available_libraries_from_parquet(parquet_dir: Path) -> Set[str]:
    """
    Get set of libraries that have permutation data available.
    
    Args:
        parquet_dir: Directory containing Parquet cluster statistics files
        
    Returns:
        Set of library names that have permutation data
    """
    if not parquet_dir.exists():
        logger.warning(f"Parquet directory not found: {parquet_dir}")
        return set()
    
    # Find any Parquet file to check library columns
    parquet_files = list(parquet_dir.glob("cluster_stats_size_*.parquet"))
    if not parquet_files:
        logger.warning(f"No Parquet files found in {parquet_dir}")
        return set()
    
    # Load one file to check available libraries
    sample_file = parquet_files[0]
    df = pd.read_parquet(sample_file)
    
    # Find boolean columns (has_*)
    bool_cols = [c for c in df.columns if c.startswith('has_')]
    
    # Also check the 'libraries' column to get actual library names
    all_libs_from_data = set()
    if 'libraries' in df.columns:
        for libs_str in df['libraries'].dropna():
            if isinstance(libs_str, str):
                libs = [lib.strip() for lib in libs_str.split(',')]
                all_libs_from_data.update(libs)
    
    # Convert column names back to library names
    available_libraries = set()
    for col in bool_cols:
        # has_go_bp -> GO BP
        # has_reactome -> Reactome
        # has_kegg -> KEGG
        lib_name = col.replace('has_', '').replace('_', ' ').title()
        # Handle special cases
        if 'go bp' in lib_name.lower():
            lib_name = 'GO BP'
        elif 'go cc' in lib_name.lower():
            lib_name = 'GO CC'
        elif 'go mf' in lib_name.lower():
            lib_name = 'GO MF'
        elif 'go ' in lib_name.lower():
            lib_name = lib_name.replace('Go ', 'GO ')
        elif 'kegg' in lib_name.lower():
            lib_name = 'KEGG'
        available_libraries.add(lib_name)
    
    # If we found libraries from the data, use those (more accurate)
    if all_libs_from_data:
        available_libraries = all_libs_from_data
    
    logger.info(f"Found {len(available_libraries)} libraries with permutation data: {sorted(available_libraries)}")
    return available_libraries


def find_intersection_libraries(
    user_selected_libraries: List[str],
    parquet_dir: Path
) -> tuple[List[str], List[str]]:
    """
    Find intersection between user-selected libraries and libraries with permutation data.
    
    Args:
        user_selected_libraries: List of library names selected by user
        parquet_dir: Directory containing Parquet cluster statistics files
        
    Returns:
        Tuple of (libraries_with_data, libraries_without_data)
    """
    available_libraries = get_available_libraries_from_parquet(parquet_dir)
    
    # Normalize library names for comparison (case-insensitive)
    available_normalized = {lib.lower().strip() for lib in available_libraries}
    user_normalized = {lib.lower().strip() for lib in user_selected_libraries}
    
    # Find intersection
    libraries_with_data = []
    libraries_without_data = []
    
    for lib in user_selected_libraries:
        lib_normalized = lib.lower().strip()
        # Try exact match first
        if lib_normalized in available_normalized:
            # Find the actual name from available libraries
            for avail_lib in available_libraries:
                if avail_lib.lower().strip() == lib_normalized:
                    libraries_with_data.append(avail_lib)
                    break
        else:
            # Try partial match (e.g., "Reactome" vs "C2: CP: Reactome Pathways")
            matched = False
            for avail_lib in available_libraries:
                if lib_normalized in avail_lib.lower() or avail_lib.lower() in lib_normalized:
                    libraries_with_data.append(avail_lib)
                    matched = True
                    break
            if not matched:
                libraries_without_data.append(lib)
    
    logger.info(f"Libraries with permutation data: {libraries_with_data}")
    if libraries_without_data:
        logger.warning(f"Libraries without permutation data (excluded from statistics): {libraries_without_data}")
    
    return libraries_with_data, libraries_without_data


def load_cluster_stats_for_size(
    parquet_dir: Path,
    gene_list_size: int,
    available_sizes: List[int]
) -> pd.DataFrame:
    """Load cluster statistics for the nearest gene list size."""
    # Round up to next increment of 50, cap at 1000
    original_size = gene_list_size
    if gene_list_size not in available_sizes:
        if gene_list_size > 1000:
            gene_list_size = 1000
            logger.warning(
                f"Gene list size {original_size} exceeds maximum permutation data size (1000). "
                f"Using size 1000 for null distribution comparison. "
                f"Note: Permutation data is only available up to size 1000."
            )
        else:
            # Round up to next multiple of 50
            gene_list_size = ((gene_list_size + 49) // 50) * 50
            logger.debug(
                f"Gene list size {original_size} not found, rounding up to nearest increment of 50: {gene_list_size}"
            )
        
        # Verify the rounded size exists
        if gene_list_size not in available_sizes:
            # Fallback: find nearest available size
            fallback_size = min(available_sizes, key=lambda x: abs(x - gene_list_size))
            logger.warning(
                f"Rounded size {gene_list_size} not found in Parquet data. "
                f"Using nearest available size: {fallback_size}"
            )
            gene_list_size = fallback_size
    
    parquet_file = parquet_dir / f"cluster_stats_size_{gene_list_size}.parquet"
    if not parquet_file.exists():
        raise FileNotFoundError(f"Parquet file not found: {parquet_file}")
    
    df = pd.read_parquet(parquet_file)
    logger.debug(f"Loaded {len(df):,} clusters from size {gene_list_size} (target: {original_size})")
    return df


def filter_by_libraries(
    df: pd.DataFrame,
    selected_libraries: List[str]
) -> pd.DataFrame:
    """
    Filter clusters by selected libraries.
    
    Includes clusters that contain at least one of the selected libraries.
    """
    if not selected_libraries:
        return df
    
    # Create filter mask: cluster must have at least one selected library
    mask = pd.Series([False] * len(df))
    
    for lib in selected_libraries:
        # Normalize library name for column matching
        col_name = f"has_{lib.lower().replace(' ', '_').replace(':', '_')}"
        if col_name in df.columns:
            mask = mask | df[col_name]
        else:
            # Fallback: check libraries string column
            mask = mask | df['libraries'].str.contains(lib, case=False, na=False)
    
    filtered_df = df[mask].copy()
    return filtered_df


def compute_null_distribution_from_raw_permutations(
    merged_permutation_file: Path,
    gene_list_size: int,
    selected_libraries: List[str],
    user_p_threshold: float,
    metrics: List[str] = None
) -> Dict[str, Dict[str, float]]:
    """
    Compute null distribution from raw permutation results, filtering by user's p-value threshold.
    
    Args:
        merged_permutation_file: Path to merged permutation results TSV file
        gene_list_size: Target gene list size
        selected_libraries: Libraries to filter by
        user_p_threshold: User's p-value threshold (must be <= 0.05)
        metrics: List of metric names to compute
        
    Returns:
        Dictionary: {gene_list_size: {metric_name: {'mean': float, 'std': float, ...}}}
    """
    from code.network_connectivity_benchmark import NetworkConnectivityAnalyzer
    
    if user_p_threshold > 0.05:
        raise ValueError(f"User p-value threshold ({user_p_threshold}) exceeds 0.05. "
                        f"Statistical benchmarking is only available for p-value thresholds <= 0.05.")
    
    if metrics is None:
        metrics = [
            'largest_cluster_genes',
            'largest_cluster_terms',
            'largest_cluster_edges',
            'largest_cluster_density',
            'largest_cluster_libraries',  # Added for library diversity benchmarking
            'largest_component_size',
            'n_connected_components',
            'avg_cluster_size',
            'avg_cluster_density',
            'avg_cluster_libraries',  # Added for library diversity benchmarking
            'fraction_in_largest_cluster',
        ]
    
    logger.info(f"Loading raw permutation results from {merged_permutation_file}")
    df = pd.read_csv(merged_permutation_file, sep='\t')
    logger.info(f"Loaded {len(df):,} permutation result rows")
    
    # Find available gene list sizes
    available_sizes = sorted(df['gene_list_size'].unique())
    logger.info(f"Available gene list sizes in permutation data: {available_sizes}")
    
    # Find matching size: round UP to next increment of 50, cap at 1000
    original_size = gene_list_size
    if gene_list_size not in available_sizes:
        # Round up to next increment of 50
        if gene_list_size > 1000:
            # Cap at 1000 and warn
            gene_list_size = 1000
            logger.warning(
                f"Gene list size {original_size} exceeds maximum permutation data size (1000). "
                f"Using size 1000 for null distribution comparison. "
                f"Note: Permutation data is only available up to size 1000."
            )
        else:
            # Round up to next multiple of 50
            gene_list_size = ((gene_list_size + 49) // 50) * 50
            size_diff = gene_list_size - original_size
            logger.info(
                f"Gene list size {original_size} not found, rounding up to nearest increment of 50: "
                f"{gene_list_size} (difference: {size_diff})"
            )
        
        # Verify the rounded size exists in available sizes
        if gene_list_size not in available_sizes:
            # Fallback: find nearest available size
            nearest_size = min(available_sizes, key=lambda x: abs(x - gene_list_size))
            logger.warning(
                f"Rounded size {gene_list_size} not found in permutation data. "
                f"Using nearest available size: {nearest_size}"
            )
            gene_list_size = nearest_size
    
    # Filter by gene list size (now using nearest)
    df = df[df['gene_list_size'] == gene_list_size].copy()
    logger.info(f"Filtered to {len(df):,} rows for gene list size {gene_list_size}")
    
    # Filter by p-value threshold
    # Handle different possible column names for p-value
    p_value_col = None
    for col in ['iteration p-value', 'p-value', 'p_value', 'P-value', 'P_value', 'Full list p-value']:
        if col in df.columns:
            p_value_col = col
            break
    
    if p_value_col is None:
        raise ValueError("Could not find p-value column in permutation results. "
                        f"Available columns: {list(df.columns)}")
    
    # Convert p-value to float and filter
    df[p_value_col] = pd.to_numeric(df[p_value_col], errors='coerce')
    df = df[df[p_value_col] <= user_p_threshold].copy()
    logger.info(f"Filtered to {len(df):,} rows with p-value <= {user_p_threshold}")
    
    # Filter by selected libraries
    if selected_libraries:
        df = df[df['Library'].isin(selected_libraries)].copy()
        logger.info(f"Filtered to {len(df):,} rows for selected libraries")
    
    # Get unique permutations
    unique_permutations = df.groupby('filename').size().reset_index()
    n_permutations = len(unique_permutations)
    logger.info(f"Processing {n_permutations} permutations")
    
    if n_permutations == 0:
        error_msg = (
            f"No permutation data available after filtering "
            f"(size={gene_list_size}, p<={user_p_threshold}, libraries={selected_libraries}). "
            f"This may occur if the p-value threshold is too strict. "
            f"Permutation data was generated with p-value threshold 0.05, so very few results may pass stricter thresholds."
        )
        raise ValueError(error_msg)
    
    # Initialize analyzer
    analyzer = NetworkConnectivityAnalyzer()
    
    # Process each permutation and compute cluster statistics
    all_cluster_stats = []
    for idx, (_, perm_info) in enumerate(unique_permutations.iterrows(), 1):
        filename = perm_info['filename']
        
        # Get data for this permutation
        perm_data = df[df['filename'] == filename]
        
        # Reset analyzer
        analyzer.reset()
        
        # Convert to iGEA results format
        results = []
        for _, row in perm_data.iterrows():
            # Try multiple possible column names for genes removed
            genes_removed = None
            for col in ['Genes removed for next iteration', 'Genes', 'genes_removed']:
                if col in row.index:
                    genes_removed = row.get(col, '')
                    break
            
            if genes_removed is None:
                genes_list = []
            elif isinstance(genes_removed, str):
                genes_list = [g.strip() for g in genes_removed.split(',') if g.strip()]
            elif isinstance(genes_removed, list):
                genes_list = genes_removed
            else:
                genes_list = []
            
            results.append({
                'Term': row.get('Term', ''),
                'Iteration': row.get('Iteration', 1),
                'Library': row.get('Library', ''),
                'Genes removed for next iteration': genes_list,
            })
        
        # Add to analyzer
        analyzer.add_igea_results(results)
        
        # Compute network-wide metrics (needed for fraction_in_largest_cluster)
        # This must be done before getting clusters to ensure metrics are available
        network_metrics = analyzer.compute_metrics(include_library_diversity=True)
        
        # Get clusters
        clusters = analyzer.get_clusters()
        
        # Store cluster statistics for this permutation
        if clusters:
            largest_cluster = clusters[0]  # Already sorted by size
            # Compute average cluster libraries
            cluster_libraries = [c.get('n_libraries', 0) for c in clusters]
            avg_cluster_libraries = np.mean(cluster_libraries) if cluster_libraries else 0.0
            
            all_cluster_stats.append({
                'filename': filename,
                'largest_cluster_genes': largest_cluster['n_genes'],
                'largest_cluster_terms': largest_cluster['n_terms'],
                'largest_cluster_edges': largest_cluster['n_edges'],
                'largest_cluster_density': largest_cluster['density'],
                'largest_cluster_libraries': largest_cluster.get('n_libraries', 0),
                'largest_component_size': largest_cluster['size'],
                'n_connected_components': len(clusters),
                'avg_cluster_size': np.mean([c['size'] for c in clusters]) if clusters else 0,
                'avg_cluster_density': np.mean([c['density'] for c in clusters]) if clusters else 0,
                'avg_cluster_libraries': avg_cluster_libraries,
                'fraction_in_largest_cluster': network_metrics.get('fraction_in_largest_cluster', 0.0),
            })
        else:
            # No clusters - all metrics are 0
            all_cluster_stats.append({
                'filename': filename,
                'largest_cluster_genes': 0,
                'largest_cluster_terms': 0,
                'largest_cluster_edges': 0,
                'largest_cluster_density': 0.0,
                'largest_cluster_libraries': 0,
                'largest_component_size': 0,
                'n_connected_components': 0,
                'avg_cluster_size': 0,
                'avg_cluster_density': 0.0,
                'avg_cluster_libraries': 0.0,
                'fraction_in_largest_cluster': 0.0,
            })
        
        if idx % 100 == 0:
            logger.info(f"Processed {idx}/{n_permutations} permutations...")
    
    # Convert to DataFrame for easier statistics computation
    stats_df = pd.DataFrame(all_cluster_stats)
    
    # Compute statistics for each metric
    null_stats = {}
    for metric in metrics:
        if metric in stats_df.columns:
            values = stats_df[metric].dropna().values
            if len(values) > 0:
                null_stats[metric] = {
                    'mean': float(np.mean(values)),
                    'std': float(np.std(values)),
                    'n': int(len(values)),
                    'min': float(np.min(values)),
                    'max': float(np.max(values)),
                    'median': float(np.median(values))
                }
            else:
                logger.warning(f"No data for metric '{metric}'")
        else:
            logger.warning(f"Metric '{metric}' not found in cluster statistics")
    
    # The actual size used (may have been rounded)
    # Use the rounded gene_list_size, not the original
    actual_size = gene_list_size  # This is the rounded size after processing
    return null_stats, actual_size  # Return both stats and the actual size used


def compute_null_distribution_from_parquet(
    parquet_dir: Path,
    gene_list_size: int,
    selected_libraries: List[str],
    metrics: List[str] = None,
    user_p_threshold: float = None,
    merged_permutation_file: Path = None
) -> Dict[str, Dict[str, float]]:
    """
    Compute null distribution statistics from Parquet cluster statistics.
    
    If user_p_threshold is provided and <= 0.05, filters raw permutation results by p-value.
    Otherwise, uses pre-computed Parquet cluster statistics (generated with p-value 0.05).
    
    Args:
        parquet_dir: Directory containing Parquet files
        gene_list_size: Target gene list size
        selected_libraries: Libraries to filter by
        metrics: List of metric names to compute (default: all cluster metrics)
        user_p_threshold: User's p-value threshold (if provided and <= 0.05, uses raw permutation filtering)
        merged_permutation_file: Path to merged permutation results TSV (required if user_p_threshold is provided)
        
    Returns:
        Dictionary: {gene_list_size: {metric_name: {'mean': float, 'std': float, ...}}}
    """
    # Check if user p-value threshold is provided and valid
    if user_p_threshold is not None:
        if user_p_threshold > 0.05:
            raise ValueError(
                f"User p-value threshold ({user_p_threshold}) exceeds 0.05. "
                f"Statistical benchmarking is only available for p-value thresholds <= 0.05. "
                f"Permutation data was generated with p-value threshold 0.05."
            )
        
        # If user threshold is exactly 0.05, use Parquet (faster, pre-computed)
        if user_p_threshold == 0.05:
            logger.info("User p-value threshold is 0.05, using pre-computed Parquet cluster statistics")
            # Continue to Parquet processing below
        else:
            # User threshold < 0.05, need to filter raw permutation results
            if merged_permutation_file is None:
                # Try to find merged permutation file
                project_root = parquet_dir.parent.parent
                possible_paths = [
                    project_root / "results" / "permutation_results" / "merged_permutation_results.tsv",
                    project_root / "archive" / "permutation_analysis" / "results" / "permutation_results" / "merged_permutation_results.tsv",
                    project_root / "merged_permutation_results.tsv",
                    parquet_dir.parent / "merged_permutation_results.tsv",
                ]
                
                merged_permutation_file = None
                for path in possible_paths:
                    if path.exists():
                        merged_permutation_file = path
                        logger.info(f"Found merged permutation file: {merged_permutation_file}")
                        break
                
                if merged_permutation_file is None:
                    raise FileNotFoundError(
                        f"Could not find merged permutation results file. "
                        f"Required for p-value filtering. Tried: {possible_paths}"
                    )
            
            logger.info(f"Using raw permutation results with p-value filtering (threshold: {user_p_threshold})")
            null_stats, actual_size = compute_null_distribution_from_raw_permutations(
                merged_permutation_file,
                gene_list_size,
                selected_libraries,
                user_p_threshold,
                metrics
            )
            # Return the actual size used (may have been rounded)
            return null_stats, actual_size
    
    # Use pre-computed Parquet cluster statistics (generated with p-value 0.05)
    if metrics is None:
        metrics = [
            'largest_cluster_genes',
            'largest_cluster_terms',
            'largest_cluster_edges',
            'largest_cluster_density',
            'largest_cluster_libraries',  # Added for library diversity benchmarking
            'largest_component_size',
            'n_connected_components',
            'avg_cluster_size',
            'avg_cluster_density',
            'avg_cluster_libraries',  # Added for library diversity benchmarking
            'fraction_in_largest_cluster',
        ]
    
    # Find available sizes
    available_sizes = sorted([
        int(f.stem.split('_')[-1])
        for f in parquet_dir.glob("cluster_stats_size_*.parquet")
    ])
    
    if not available_sizes:
        raise FileNotFoundError(f"No Parquet files found in {parquet_dir}")
    
    # Load cluster statistics (this will round up to nearest increment of 50)
    original_size = gene_list_size
    df = load_cluster_stats_for_size(parquet_dir, gene_list_size, available_sizes)
    
    # Get the actual size used (may have been rounded)
    # The load_cluster_stats_for_size function rounds internally, so we need to determine
    # what size was actually loaded by checking what parquet file exists
    actual_size = original_size
    if original_size not in available_sizes:
        if original_size > 1000:
            actual_size = 1000
        else:
            actual_size = ((original_size + 49) // 50) * 50
        if actual_size not in available_sizes:
            actual_size = min(available_sizes, key=lambda x: abs(x - actual_size))
        if actual_size != original_size:
            logger.info(f"Using null distribution from size {actual_size} (rounded from {original_size})")
    
    # Filter by libraries
    if selected_libraries:
        df = filter_by_libraries(df, selected_libraries)
        logger.info(f"Filtered to {len(df):,} clusters containing selected libraries")
    
    # Compute statistics for each metric
    # All metrics need to be aggregated per permutation (one value per permutation)
    null_stats = {}
    
    # Helper function to get largest cluster metric per permutation
    def get_largest_cluster_metric(group, metric_col):
        """Get metric from largest cluster (cluster_number=1) for a permutation."""
        if len(group) == 0:
            return 0 if metric_col != 'density' else 0.0
        # cluster_number=1 is the largest cluster
        largest_clusters = group[group['cluster_number'] == 1]
        if len(largest_clusters) > 0:
            return largest_clusters.iloc[0][metric_col]
        # Fallback: get first cluster if no cluster_number=1
        return group.iloc[0][metric_col]
    
    for metric in metrics:
        values = None
        
        if metric == 'largest_cluster_genes':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'n_genes')
            )
            values = largest_per_perm.values
        elif metric == 'largest_cluster_terms':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'n_terms')
            )
            values = largest_per_perm.values
        elif metric == 'largest_cluster_edges':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'n_edges')
            )
            values = largest_per_perm.values
        elif metric == 'largest_cluster_density':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'density')
            )
            values = largest_per_perm.values
        elif metric == 'largest_component_size':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'cluster_size')
            )
            values = largest_per_perm.values
        elif metric == 'n_connected_components':
            # Count clusters per permutation
            clusters_per_perm = df.groupby('filename').size()
            values = clusters_per_perm.values
        elif metric == 'avg_cluster_size':
            avg_per_perm = df.groupby('filename')['cluster_size'].mean()
            values = avg_per_perm.values
        elif metric == 'avg_cluster_density':
            avg_per_perm = df.groupby('filename')['density'].mean()
            values = avg_per_perm.values
        elif metric == 'weighted_avg_cluster_density':
            # Weighted by cluster size
            def weighted_density(group):
                if len(group) == 0:
                    return 0.0
                total_size = group['cluster_size'].sum()
                if total_size == 0:
                    return 0.0
                return (group['density'] * group['cluster_size']).sum() / total_size
            weighted_per_perm = df.groupby('filename', group_keys=False).apply(weighted_density)
            values = weighted_per_perm.values
        elif metric == 'fraction_in_largest_cluster':
            def compute_fraction(group):
                if len(group) == 0:
                    return 0.0
                largest = group[group['cluster_number'] == 1]
                if len(largest) == 0:
                    largest = group.iloc[[0]]
                total_genes = group['n_genes'].sum()
                total_terms = group['n_terms'].sum()
                total_nodes = total_genes + total_terms
                if total_nodes == 0:
                    return 0.0
                return largest.iloc[0]['cluster_size'] / total_nodes
            fractions = df.groupby('filename', group_keys=False).apply(compute_fraction)
            values = fractions.values
        elif metric == 'fraction_edges_in_largest_cluster':
            def compute_edge_fraction(group):
                if len(group) == 0:
                    return 0.0
                largest = group[group['cluster_number'] == 1]
                if len(largest) == 0:
                    largest = group.iloc[[0]]
                total_edges = group['n_edges'].sum()
                if total_edges == 0:
                    return 0.0
                return largest.iloc[0]['n_edges'] / total_edges
            edge_fractions = df.groupby('filename', group_keys=False).apply(compute_edge_fraction)
            values = edge_fractions.values
        elif metric == 'largest_cluster_libraries':
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, 'n_libraries')
            )
            values = largest_per_perm.values
        elif metric == 'avg_cluster_libraries':
            avg_per_perm = df.groupby('filename')['n_libraries'].mean()
            values = avg_per_perm.values
        elif metric in ['cluster_size', 'n_genes', 'n_terms', 'n_edges', 'density', 'n_libraries']:
            # These are cluster-level metrics - aggregate per permutation (use largest cluster)
            largest_per_perm = df.groupby('filename', group_keys=False).apply(
                lambda x: get_largest_cluster_metric(x, metric)
            )
            values = largest_per_perm.values
        else:
            logger.warning(f"Metric '{metric}' not found and no computation method available")
            continue
        
        if values is not None and len(values) > 0:
            null_stats[metric] = {
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'n': int(len(values)),
                'min': float(np.min(values)),
                'max': float(np.max(values)),
                'median': float(np.median(values))
            }
        else:
            logger.warning(f"No data for metric '{metric}'")
    
    return null_stats, actual_size  # Return both stats and the actual size used


def compute_null_distribution_parallel(
    parquet_dir: Path,
    gene_list_size: int,
    selected_libraries: List[str],
    result_dict: Dict,
    lock: threading.Lock,
    user_p_threshold: float = None,
    merged_permutation_file: Path = None
) -> None:
    """
    Compute null distribution in a separate thread and store result in shared dictionary.
    
    Args:
        parquet_dir: Directory containing Parquet files
        gene_list_size: Target gene list size
        selected_libraries: Libraries to filter by
        result_dict: Shared dictionary to store results
        lock: Thread lock for safe dictionary access
        user_p_threshold: User's p-value threshold (if > 0.05, benchmarking unavailable)
        merged_permutation_file: Path to merged permutation results TSV (for p-value filtering)
    """
    try:
        logger.info(f"Computing null distribution for size {gene_list_size} with libraries: {selected_libraries}")
        
        # Check if user p-value threshold is too high
        if user_p_threshold is not None and user_p_threshold > 0.05:
            with lock:
                result_dict['status'] = 'unavailable'
                result_dict['error'] = (
                    f"Statistical benchmarking is not available for p-value thresholds > 0.05. "
                    f"Your threshold: {user_p_threshold}. "
                    f"Permutation data was generated with p-value threshold 0.05."
                )
            logger.warning(result_dict['error'])
            return
        
        if user_p_threshold is not None:
            logger.info(f"User p-value threshold: {user_p_threshold} (filtering permutation data)")
        else:
            logger.info(f"Using pre-computed cluster statistics (permutation data generated with 0.05)")
        
        # Compute null distribution (may round up the size)
        original_size = gene_list_size
        null_stats, actual_size = compute_null_distribution_from_parquet(
            parquet_dir,
            gene_list_size,
            selected_libraries,
            user_p_threshold=user_p_threshold,
            merged_permutation_file=merged_permutation_file
        )
        
        with lock:
            # Store with the actual size used (may be rounded up)
            result_dict['null_distribution'] = {str(actual_size): null_stats}
            result_dict['original_gene_list_size'] = original_size
            result_dict['null_distribution_size'] = actual_size
            result_dict['status'] = 'completed'
            result_dict['libraries_used'] = selected_libraries
            result_dict['permutation_p_threshold'] = 0.05  # Permutations always generated with 0.05
            if user_p_threshold is not None:
                result_dict['user_p_threshold'] = user_p_threshold
        
        logger.info(f"✓ Null distribution computation completed")
    except ValueError as e:
        # Handle p-value threshold errors
        logger.error(f"Error: {e}")
        with lock:
            result_dict['status'] = 'unavailable'
            result_dict['error'] = str(e)
    except Exception as e:
        logger.error(f"Error computing null distribution: {e}", exc_info=True)
        with lock:
            result_dict['status'] = 'error'
            result_dict['error'] = str(e)
