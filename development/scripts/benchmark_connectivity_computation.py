#!/usr/bin/env python3
"""
Benchmark connectivity computation time for a gene list.

This script:
1. Loads pre-segmented Parquet cluster statistics
2. Filters by selected libraries
3. Computes null distribution statistics
4. Times the operation (without caching)

Usage:
    python3 benchmark_connectivity_computation.py --gene-list-size 200 --libraries Reactome "GO BP" "GO CC" "GO MF" KEGG
"""

import pandas as pd
import numpy as np
import time
import logging
from pathlib import Path
from typing import List, Dict, Optional
import json

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def find_nearest_gene_list_sizes(target_size: int, available_sizes: List[int]) -> List[int]:
    """
    Find the nearest gene list size(s) in available data.
    
    Returns list of sizes (usually 1, sometimes 2 for interpolation).
    """
    available_sizes = sorted(available_sizes)
    
    # Find nearest size
    nearest = min(available_sizes, key=lambda x: abs(x - target_size))
    
    # If target is between two sizes, return both for potential interpolation
    if target_size < nearest:
        smaller = max([s for s in available_sizes if s < target_size], default=None)
        if smaller:
            return [smaller, nearest]
    elif target_size > nearest:
        larger = min([s for s in available_sizes if s > target_size], default=None)
        if larger:
            return [nearest, larger]
    
    return [nearest]


def load_cluster_stats_for_size(
    parquet_dir: Path,
    gene_list_size: int,
    available_sizes: List[int]
) -> pd.DataFrame:
    """Load cluster statistics for the nearest gene list size."""
    nearest_sizes = find_nearest_gene_list_sizes(gene_list_size, available_sizes)
    
    logger.info(f"Target size: {gene_list_size}, Loading sizes: {nearest_sizes}")
    
    dfs = []
    for size in nearest_sizes:
        parquet_file = parquet_dir / f"cluster_stats_size_{size}.parquet"
        if parquet_file.exists():
            df = pd.read_parquet(parquet_file)
            logger.info(f"  Loaded size {size}: {len(df):,} clusters")
            dfs.append(df)
        else:
            logger.warning(f"  File not found: {parquet_file}")
    
    if not dfs:
        raise FileNotFoundError(f"No cluster statistics found for size {gene_list_size}")
    
    # Combine if multiple sizes
    if len(dfs) > 1:
        combined_df = pd.concat(dfs, ignore_index=True)
        logger.info(f"  Combined: {len(combined_df):,} clusters")
        return combined_df
    else:
        return dfs[0]


def filter_by_libraries(
    df: pd.DataFrame,
    selected_libraries: List[str]
) -> pd.DataFrame:
    """
    Filter clusters by selected libraries.
    
    Includes clusters that contain at least one of the selected libraries.
    """
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


def compute_null_distribution_statistics(
    df: pd.DataFrame,
    metrics: List[str]
) -> Dict[str, Dict[str, float]]:
    """
    Compute mean and std for each metric from cluster statistics.
    
    Args:
        df: Filtered cluster statistics DataFrame
        metrics: List of metric column names
        
    Returns:
        Dictionary: {metric_name: {'mean': float, 'std': float, 'n': int}}
    """
    null_stats = {}
    
    for metric in metrics:
        if metric not in df.columns:
            logger.warning(f"Metric '{metric}' not found in DataFrame, skipping")
            continue
        
        values = df[metric].dropna()
        if len(values) > 0:
            null_stats[metric] = {
                'mean': float(values.mean()),
                'std': float(values.std()),
                'n': int(len(values)),
                'min': float(values.min()),
                'max': float(values.max()),
                'median': float(values.median())
            }
        else:
            logger.warning(f"No data for metric '{metric}'")
    
    return null_stats


def benchmark_connectivity_computation(
    parquet_dir: str,
    gene_list_size: int,
    selected_libraries: List[str],
    metrics: List[str] = None
) -> Dict:
    """
    Benchmark the connectivity computation process.
    
    Returns:
        Dictionary with timing information and computed statistics
    """
    parquet_path = Path(parquet_dir)
    
    if not parquet_path.exists():
        raise FileNotFoundError(f"Parquet directory not found: {parquet_dir}")
    
    # Default metrics to compute
    if metrics is None:
        metrics = [
            'cluster_size',
            'n_genes',
            'n_terms',
            'n_edges',
            'density',
            'n_libraries'
        ]
    
    logger.info("=" * 80)
    logger.info("BENCHMARK: Connectivity Computation")
    logger.info("=" * 80)
    logger.info(f"Gene list size: {gene_list_size}")
    logger.info(f"Selected libraries: {selected_libraries}")
    logger.info(f"Metrics: {metrics}")
    logger.info("=" * 80)
    
    # Find available sizes
    available_sizes = sorted([
        int(f.stem.split('_')[-1])
        for f in parquet_path.glob("cluster_stats_size_*.parquet")
    ])
    logger.info(f"Available gene list sizes: {available_sizes}")
    
    # Time the entire process
    start_time = time.time()
    
    # Step 1: Load Parquet file(s)
    load_start = time.time()
    df = load_cluster_stats_for_size(parquet_path, gene_list_size, available_sizes)
    load_time = time.time() - load_start
    logger.info(f"✓ Load time: {load_time:.3f} seconds ({len(df):,} clusters)")
    
    # Step 2: Filter by libraries
    filter_start = time.time()
    filtered_df = filter_by_libraries(df, selected_libraries)
    filter_time = time.time() - filter_start
    logger.info(f"✓ Filter time: {filter_time:.3f} seconds ({len(filtered_df):,} clusters after filtering)")
    
    # Step 3: Compute statistics
    stats_start = time.time()
    null_stats = compute_null_distribution_statistics(filtered_df, metrics)
    stats_time = time.time() - stats_start
    logger.info(f"✓ Statistics computation time: {stats_time:.3f} seconds")
    
    total_time = time.time() - start_time
    
    logger.info("=" * 80)
    logger.info(f"TOTAL TIME: {total_time:.3f} seconds")
    logger.info("=" * 80)
    logger.info("\nBreakdown:")
    logger.info(f"  Load:        {load_time:.3f}s ({load_time/total_time*100:.1f}%)")
    logger.info(f"  Filter:      {filter_time:.3f}s ({filter_time/total_time*100:.1f}%)")
    logger.info(f"  Statistics:  {stats_time:.3f}s ({stats_time/total_time*100:.1f}%)")
    logger.info("=" * 80)
    
    # Show sample statistics
    logger.info("\nSample null distribution statistics:")
    for metric, stats in list(null_stats.items())[:3]:
        logger.info(f"  {metric}: mean={stats['mean']:.2f}, std={stats['std']:.2f}, n={stats['n']}")
    
    return {
        'total_time': total_time,
        'load_time': load_time,
        'filter_time': filter_time,
        'stats_time': stats_time,
        'n_clusters_loaded': len(df),
        'n_clusters_filtered': len(filtered_df),
        'null_statistics': null_stats,
        'gene_list_size': gene_list_size,
        'selected_libraries': selected_libraries
    }


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Benchmark connectivity computation time'
    )
    parser.add_argument(
        '--parquet-dir',
        type=str,
        default='results/permutation_cluster_statistics_parquet',
        help='Directory containing Parquet cluster statistics files'
    )
    parser.add_argument(
        '--gene-list-size',
        type=int,
        default=200,
        help='Gene list size to benchmark'
    )
    parser.add_argument(
        '--libraries',
        type=str,
        nargs='+',
        default=['Reactome', 'KEGG', 'GO BP', 'GO MF', 'GO CC'],
        help='Libraries to filter by'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Optional: Save benchmark results to JSON file'
    )
    
    args = parser.parse_args()
    
    result = benchmark_connectivity_computation(
        args.parquet_dir,
        args.gene_list_size,
        args.libraries
    )
    
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        logger.info(f"\n✓ Results saved to {args.output}")


if __name__ == '__main__':
    main()
