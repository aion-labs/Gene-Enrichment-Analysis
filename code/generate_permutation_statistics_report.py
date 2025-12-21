#!/usr/bin/env python3
"""
Generate statistical reports for permutation data across gene list sizes and p-value thresholds.

This script:
1. Loads merged permutation results
2. For each gene list size (50-1000) and p-value cutoff (0.05, 0.01, 0.001):
   - Filters permutation data
   - Computes clusters for each permutation
   - Extracts largest cluster metrics
   - Computes mean and standard deviation for each metric
3. Generates a summary report (text file)
4. Generates a detailed CSV with all clusters (one row per cluster per permutation)
"""

import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np
from collections import defaultdict

# Add code directory to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from network_connectivity_benchmark import NetworkConnectivityAnalyzer

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# Gene list sizes to process
GENE_LIST_SIZES = [50, 100, 150, 200, 250, 300, 400, 500, 750, 1000]

# P-value thresholds to process
P_VALUE_THRESHOLDS = [0.05, 0.01, 0.001]

# Metrics to compute for largest cluster
LARGEST_CLUSTER_METRICS = [
    'largest_cluster_genes',
    'largest_cluster_terms',
    'largest_cluster_edges',
    'largest_cluster_density',
    'largest_component_size',  # Total size (genes + terms)
    'fraction_in_largest_cluster',
    'fraction_edges_in_largest_cluster',
    'largest_cluster_libraries',
]


def find_merged_permutation_file() -> Path:
    """Find the merged permutation results file."""
    possible_paths = [
        PROJECT_ROOT / "archive" / "permutation_analysis" / "results" / "permutation_results" / "merged_permutation_results.tsv",
        PROJECT_ROOT / "results" / "permutation_results" / "merged_permutation_results.tsv",
        PROJECT_ROOT / "merged_permutation_results.tsv",
    ]
    
    for path in possible_paths:
        if path.exists():
            logger.info(f"Found merged permutation file: {path}")
            return path
    
    raise FileNotFoundError(
        f"Could not find merged permutation results file. Tried:\n" +
        "\n".join(f"  - {p}" for p in possible_paths)
    )


def find_p_value_column(df: pd.DataFrame) -> str:
    """Find the p-value column name in the dataframe."""
    possible_names = [
        'iteration p-value', 'p-value', 'p_value', 'P-value', 'P_value', 
        'Full list p-value', 'iteration_p-value'
    ]
    
    for col in possible_names:
        if col in df.columns:
            return col
    
    raise ValueError(
        f"Could not find p-value column. Available columns: {list(df.columns)}"
    )


def process_permutation_data(
    df: pd.DataFrame,
    gene_list_size: int,
    p_value_threshold: float,
    p_value_col: str
) -> Tuple[List[Dict], int]:
    """
    Process permutation data for a specific size and p-value threshold.
    
    Returns:
        (list of cluster records, number of permutations processed)
    """
    # Filter by gene list size
    size_df = df[df['gene_list_size'] == gene_list_size].copy()
    if len(size_df) == 0:
        logger.warning(f"No data for gene list size {gene_list_size}")
        return [], 0
    
    # Filter by p-value threshold
    size_df[p_value_col] = pd.to_numeric(size_df[p_value_col], errors='coerce')
    filtered_df = size_df[size_df[p_value_col] <= p_value_threshold].copy()
    
    if len(filtered_df) == 0:
        logger.warning(
            f"No data for gene list size {gene_list_size} with p-value <= {p_value_threshold}"
        )
        return [], 0
    
    # Get unique permutations
    unique_permutations = filtered_df['filename'].unique()
    n_permutations = len(unique_permutations)
    
    logger.info(
        f"Processing {n_permutations} permutations for size {gene_list_size}, "
        f"p-value <= {p_value_threshold}"
    )
    
    # Initialize analyzer
    analyzer = NetworkConnectivityAnalyzer()
    
    # Process each permutation
    all_cluster_records = []
    
    for perm_idx, filename in enumerate(unique_permutations, 1):
        # Get data for this permutation
        perm_data = filtered_df[filtered_df['filename'] == filename].copy()
        
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
        
        # Get clusters
        clusters = analyzer.get_clusters()
        
        # Compute network-wide metrics once per permutation (needed for largest cluster)
        metrics = None
        if clusters:
            metrics = analyzer.compute_metrics(include_library_diversity=True)
        
        # Extract permutation ID from filename
        # Format might be: permutation_001.tsv, random_200_genes_001.tsv, etc.
        perm_id = Path(filename).stem if filename else f"permutation_{perm_idx}"
        
        # Store all clusters for this permutation
        if clusters:
            for cluster in clusters:
                cluster_record = {
                    'gene_list_size': gene_list_size,
                    'p_value_threshold': p_value_threshold,
                    'permutation_id': perm_id,
                    'filename': filename,
                    'cluster_number': cluster.get('cluster_number', 0),
                    'n_genes': cluster['n_genes'],
                    'n_terms': cluster['n_terms'],
                    'n_edges': cluster['n_edges'],
                    'cluster_size': cluster['size'],  # Total size (genes + terms)
                    'density': cluster['density'],
                    'n_libraries': cluster.get('n_libraries', 0),
                }
                
                # Add network-wide metrics for largest cluster only
                if cluster.get('cluster_number', 0) == 1 and metrics:  # Largest cluster
                    cluster_record.update({
                        'largest_cluster_genes': cluster['n_genes'],
                        'largest_cluster_terms': cluster['n_terms'],
                        'largest_cluster_edges': cluster['n_edges'],
                        'largest_cluster_density': cluster['density'],
                        'largest_component_size': cluster['size'],
                        'fraction_in_largest_cluster': metrics.get('fraction_in_largest_cluster', 0.0),
                        'fraction_edges_in_largest_cluster': metrics.get('fraction_edges_in_largest_cluster', 0.0),
                        'largest_cluster_libraries': cluster.get('n_libraries', 0),
                        'n_connected_components': metrics.get('n_connected_components', 0),
                        'n_genes_total': metrics.get('n_genes', 0),
                        'n_terms_total': metrics.get('n_terms', 0),
                        'n_edges_total': metrics.get('n_edges', 0),
                    })
                else:
                    # For non-largest clusters, set these to None/0
                    cluster_record.update({
                        'largest_cluster_genes': None,
                        'largest_cluster_terms': None,
                        'largest_cluster_edges': None,
                        'largest_cluster_density': None,
                        'largest_component_size': None,
                        'fraction_in_largest_cluster': None,
                        'fraction_edges_in_largest_cluster': None,
                        'largest_cluster_libraries': None,
                        'n_connected_components': None,
                        'n_genes_total': None,
                        'n_terms_total': None,
                        'n_edges_total': None,
                    })
                
                all_cluster_records.append(cluster_record)
        else:
            # No clusters - create a record with zeros
            cluster_record = {
                'gene_list_size': gene_list_size,
                'p_value_threshold': p_value_threshold,
                'permutation_id': perm_id,
                'filename': filename,
                'cluster_number': 0,
                'n_genes': 0,
                'n_terms': 0,
                'n_edges': 0,
                'cluster_size': 0,
                'density': 0.0,
                'n_libraries': 0,
                'largest_cluster_genes': 0,
                'largest_cluster_terms': 0,
                'largest_cluster_edges': 0,
                'largest_cluster_density': 0.0,
                'largest_component_size': 0,
                'fraction_in_largest_cluster': 0.0,
                'fraction_edges_in_largest_cluster': 0.0,
                'largest_cluster_libraries': 0,
                'n_connected_components': 0,
                'n_genes_total': 0,
                'n_terms_total': 0,
                'n_edges_total': 0,
            }
            all_cluster_records.append(cluster_record)
        
        if perm_idx % 100 == 0:
            logger.info(f"  Processed {perm_idx}/{n_permutations} permutations...")
    
    return all_cluster_records, n_permutations


def compute_statistics_for_largest_clusters(
    cluster_df: pd.DataFrame,
    gene_list_size: int,
    p_value_threshold: float
) -> Dict[str, Dict[str, float]]:
    """
    Compute statistics (mean, std) for largest cluster metrics.
    
    Only includes the largest cluster (cluster_number == 1) from each permutation.
    """
    # Filter to largest clusters only
    largest_clusters = cluster_df[
        (cluster_df['gene_list_size'] == gene_list_size) &
        (cluster_df['p_value_threshold'] == p_value_threshold) &
        (cluster_df['cluster_number'] == 1)
    ].copy()
    
    if len(largest_clusters) == 0:
        return {}
    
    # Compute statistics for each metric
    stats = {}
    
    for metric in LARGEST_CLUSTER_METRICS:
        if metric in largest_clusters.columns:
            values = largest_clusters[metric].dropna()
            if len(values) > 0:
                stats[metric] = {
                    'mean': float(np.mean(values)),
                    'std': float(np.std(values)),
                    'n': int(len(values)),
                    'min': float(np.min(values)),
                    'max': float(np.max(values)),
                    'median': float(np.median(values)),
                }
            else:
                stats[metric] = {
                    'mean': 0.0,
                    'std': 0.0,
                    'n': 0,
                    'min': 0.0,
                    'max': 0.0,
                    'median': 0.0,
                }
    
    return stats


def generate_summary_report(
    all_stats: Dict[Tuple[int, float], Dict[str, Dict[str, float]]],
    output_file: Path
) -> None:
    """Generate a text summary report."""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Permutation Statistics Report\n")
        f.write("=" * 80 + "\n\n")
        f.write("This report contains mean and standard deviation for largest cluster metrics\n")
        f.write("across different gene list sizes and p-value thresholds.\n\n")
        f.write("=" * 80 + "\n\n")
        
        # Group by gene list size
        for gene_list_size in sorted(GENE_LIST_SIZES):
            f.write(f"\n{'=' * 80}\n")
            f.write(f"Gene List Size: {gene_list_size}\n")
            f.write(f"{'=' * 80}\n\n")
            
            for p_value in sorted(P_VALUE_THRESHOLDS, reverse=True):
                key = (gene_list_size, p_value)
                if key not in all_stats or not all_stats[key]:
                    f.write(f"P-value threshold: {p_value}\n")
                    f.write(f"  No data available\n\n")
                    continue
                
                stats = all_stats[key]
                n_permutations = stats.get('largest_cluster_genes', {}).get('n', 0)
                
                f.write(f"P-value threshold: {p_value}\n")
                f.write(f"  Number of permutations: {n_permutations}\n\n")
                
                # Write statistics for each metric
                for metric in LARGEST_CLUSTER_METRICS:
                    if metric in stats:
                        metric_stats = stats[metric]
                        f.write(f"  {metric}:\n")
                        f.write(f"    Mean:   {metric_stats['mean']:.4f}\n")
                        f.write(f"    Std:    {metric_stats['std']:.4f}\n")
                        f.write(f"    Min:    {metric_stats['min']:.4f}\n")
                        f.write(f"    Max:    {metric_stats['max']:.4f}\n")
                        f.write(f"    Median: {metric_stats['median']:.4f}\n")
                        f.write(f"    N:      {metric_stats['n']}\n")
                        f.write("\n")
                
                f.write("\n")
    
    logger.info(f"Summary report written to: {output_file}")


def generate_summary_csv(
    all_stats: Dict[Tuple[int, float], Dict[str, Dict[str, float]]],
    output_file: Path
) -> None:
    """Generate a CSV summary with mean and std for each size/p-value combination."""
    rows = []
    
    for gene_list_size in sorted(GENE_LIST_SIZES):
        for p_value in sorted(P_VALUE_THRESHOLDS, reverse=True):
            key = (gene_list_size, p_value)
            if key not in all_stats or not all_stats[key]:
                continue
            
            stats = all_stats[key]
            n_permutations = stats.get('largest_cluster_genes', {}).get('n', 0)
            
            row = {
                'gene_list_size': gene_list_size,
                'p_value_threshold': p_value,
                'n_permutations': n_permutations,
            }
            
            # Add mean and std for each metric
            for metric in LARGEST_CLUSTER_METRICS:
                if metric in stats:
                    row[f'{metric}_mean'] = stats[metric]['mean']
                    row[f'{metric}_std'] = stats[metric]['std']
                    row[f'{metric}_min'] = stats[metric]['min']
                    row[f'{metric}_max'] = stats[metric]['max']
                    row[f'{metric}_median'] = stats[metric]['median']
                else:
                    row[f'{metric}_mean'] = None
                    row[f'{metric}_std'] = None
                    row[f'{metric}_min'] = None
                    row[f'{metric}_max'] = None
                    row[f'{metric}_median'] = None
            
            rows.append(row)
    
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(output_file, index=False)
    logger.info(f"Summary CSV written to: {output_file}")


def main():
    """Main function."""
    # Find merged permutation file
    merged_file = find_merged_permutation_file()
    
    # Load permutation data
    logger.info(f"Loading permutation data from {merged_file}...")
    df = pd.read_csv(merged_file, sep='\t')
    logger.info(f"Loaded {len(df):,} rows")
    
    # Find p-value column
    p_value_col = find_p_value_column(df)
    logger.info(f"Using p-value column: {p_value_col}")
    
    # Check available gene list sizes
    available_sizes = sorted(df['gene_list_size'].unique())
    logger.info(f"Available gene list sizes: {available_sizes}")
    
    # Process all combinations
    all_cluster_records = []
    all_stats = {}
    
    for gene_list_size in GENE_LIST_SIZES:
        if gene_list_size not in available_sizes:
            logger.warning(f"Skipping gene list size {gene_list_size} (not in data)")
            continue
        
        for p_value in P_VALUE_THRESHOLDS:
            logger.info(f"\n{'=' * 80}")
            logger.info(f"Processing: Size {gene_list_size}, P-value <= {p_value}")
            logger.info(f"{'=' * 80}")
            
            # Process permutation data
            cluster_records, n_permutations = process_permutation_data(
                df, gene_list_size, p_value, p_value_col
            )
            
            if n_permutations == 0:
                logger.warning(f"No permutations found for size {gene_list_size}, p-value <= {p_value}")
                continue
            
            all_cluster_records.extend(cluster_records)
    
    # Create DataFrame from all cluster records
    logger.info(f"\nCreating cluster DataFrame...")
    cluster_df = pd.DataFrame(all_cluster_records)
    logger.info(f"Total cluster records: {len(cluster_df):,}")
    
    # Compute statistics for each size/p-value combination
    logger.info(f"\nComputing statistics...")
    for gene_list_size in GENE_LIST_SIZES:
        for p_value in P_VALUE_THRESHOLDS:
            stats = compute_statistics_for_largest_clusters(
                cluster_df, gene_list_size, p_value
            )
            if stats:
                all_stats[(gene_list_size, p_value)] = stats
    
    # Create output directory
    output_dir = PROJECT_ROOT / "results" / "permutation_statistics"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate summary report
    summary_report = output_dir / "permutation_statistics_report.txt"
    generate_summary_report(all_stats, summary_report)
    
    # Generate summary CSV
    summary_csv = output_dir / "permutation_statistics_summary.csv"
    generate_summary_csv(all_stats, summary_csv)
    
    # Generate detailed cluster CSV
    detailed_csv = output_dir / "permutation_clusters_detailed.csv"
    cluster_df.to_csv(detailed_csv, index=False)
    logger.info(f"Detailed cluster CSV written to: {detailed_csv}")
    logger.info(f"  Total clusters: {len(cluster_df):,}")
    logger.info(f"  Columns: {list(cluster_df.columns)}")
    
    logger.info(f"\n{'=' * 80}")
    logger.info("Processing complete!")
    logger.info(f"{'=' * 80}")
    logger.info(f"Output files:")
    logger.info(f"  1. Summary report: {summary_report}")
    logger.info(f"  2. Summary CSV: {summary_csv}")
    logger.info(f"  3. Detailed clusters CSV: {detailed_csv}")


if __name__ == "__main__":
    main()

