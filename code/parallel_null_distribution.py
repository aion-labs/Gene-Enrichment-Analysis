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
    # Find nearest size
    nearest = min(available_sizes, key=lambda x: abs(x - gene_list_size))
    
    parquet_file = parquet_dir / f"cluster_stats_size_{nearest}.parquet"
    if not parquet_file.exists():
        raise FileNotFoundError(f"Parquet file not found: {parquet_file}")
    
    df = pd.read_parquet(parquet_file)
    logger.debug(f"Loaded {len(df):,} clusters from size {nearest} (target: {gene_list_size})")
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


def compute_null_distribution_from_parquet(
    parquet_dir: Path,
    gene_list_size: int,
    selected_libraries: List[str],
    metrics: List[str] = None
) -> Dict[str, Dict[str, float]]:
    """
    Compute null distribution statistics from Parquet cluster statistics.
    
    Args:
        parquet_dir: Directory containing Parquet files
        gene_list_size: Target gene list size
        selected_libraries: Libraries to filter by
        metrics: List of metric names to compute (default: all cluster metrics)
        
    Returns:
        Dictionary: {gene_list_size: {metric_name: {'mean': float, 'std': float, ...}}}
    """
    if metrics is None:
        metrics = [
            'cluster_size',
            'n_genes',
            'n_terms',
            'n_edges',
            'density',
            'n_libraries',
            'largest_cluster_genes',
            'largest_cluster_terms',
            'largest_cluster_edges',
            'largest_cluster_density',
            'largest_component_size',
            'n_connected_components',
            'avg_cluster_size',
            'avg_cluster_density',
            'fraction_in_largest_cluster',
        ]
    
    # Find available sizes
    available_sizes = sorted([
        int(f.stem.split('_')[-1])
        for f in parquet_dir.glob("cluster_stats_size_*.parquet")
    ])
    
    if not available_sizes:
        raise FileNotFoundError(f"No Parquet files found in {parquet_dir}")
    
    # Load cluster statistics
    df = load_cluster_stats_for_size(parquet_dir, gene_list_size, available_sizes)
    
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
    
    return null_stats


def compute_null_distribution_parallel(
    parquet_dir: Path,
    gene_list_size: int,
    selected_libraries: List[str],
    result_dict: Dict,
    lock: threading.Lock
) -> None:
    """
    Compute null distribution in a separate thread and store result in shared dictionary.
    
    Args:
        parquet_dir: Directory containing Parquet files
        gene_list_size: Target gene list size
        selected_libraries: Libraries to filter by
        result_dict: Shared dictionary to store results
        lock: Thread lock for safe dictionary access
    """
    try:
        logger.info(f"Computing null distribution for size {gene_list_size} with libraries: {selected_libraries}")
        null_stats = compute_null_distribution_from_parquet(
            parquet_dir,
            gene_list_size,
            selected_libraries
        )
        
        with lock:
            result_dict['null_distribution'] = {str(gene_list_size): null_stats}
            result_dict['status'] = 'completed'
            result_dict['libraries_used'] = selected_libraries
        
        logger.info(f"✓ Null distribution computation completed")
    except Exception as e:
        logger.error(f"Error computing null distribution: {e}", exc_info=True)
        with lock:
            result_dict['status'] = 'error'
            result_dict['error'] = str(e)
