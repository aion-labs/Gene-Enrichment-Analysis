#!/usr/bin/env python3
"""
Test script to compute clusters from random permutations with p-value filtering.

This script:
1. Loads raw permutation results
2. Filters by p-value threshold
3. Computes clusters for each permutation
4. Shows statistics
"""

import sys
import logging
from pathlib import Path
import pandas as pd
from typing import Dict, List

# Add code directory to path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from network_connectivity_benchmark import NetworkConnectivityAnalyzer
from parallel_null_distribution import compute_null_distribution_from_raw_permutations

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Paths
RESULTS_DIR = PROJECT_ROOT / "results"
PARQUET_DIR = RESULTS_DIR / "permutation_cluster_statistics_parquet"
# Try multiple possible locations for merged permutation file
MERGED_PERMUTATION_FILE = None
possible_paths = [
    RESULTS_DIR / "permutation_results" / "merged_permutation_results.tsv",
    PROJECT_ROOT / "archive" / "permutation_analysis" / "results" / "permutation_results" / "merged_permutation_results.tsv",
    PROJECT_ROOT / "merged_permutation_results.tsv",
]
for path in possible_paths:
    if path.exists():
        MERGED_PERMUTATION_FILE = path
        break

# Test parameters
TEST_GENE_LIST_SIZE = 200
TEST_P_THRESHOLD = 0.05
TEST_LIBRARIES = ['Reactome', 'KEGG', 'GO BP', 'GO MF', 'GO CC']


def test_cluster_computation_from_permutations():
    """Test computing clusters from filtered permutation data."""
    print("\n" + "=" * 80)
    print("TEST: Cluster Computation from Random Permutations")
    print("=" * 80)
    
    # Check if merged permutation file exists
    if MERGED_PERMUTATION_FILE is None or not MERGED_PERMUTATION_FILE.exists():
        print(f"\n❌ Merged permutation file not found")
        print(f"   Tried paths:")
        for path in possible_paths:
            print(f"     - {path} ({'exists' if path.exists() else 'not found'})")
        print("   Please ensure permutation results are available.")
        return
    
    print(f"\n[1/4] Loading permutation results...")
    print(f"   File: {MERGED_PERMUTATION_FILE}")
    df = pd.read_csv(MERGED_PERMUTATION_FILE, sep='\t')
    print(f"   ✓ Loaded {len(df):,} permutation result rows")
    print(f"   Columns: {list(df.columns)}")
    
    # Check for p-value column (try multiple possible names)
    p_value_col = None
    for col in ['iteration p-value', 'p-value', 'p_value', 'P-value', 'P_value', 'Full list p-value']:
        if col in df.columns:
            p_value_col = col
            break
    
    if p_value_col is None:
        print(f"\n❌ Could not find p-value column in permutation results")
        print(f"   Available columns: {list(df.columns)}")
        return
    
    print(f"   ✓ Found p-value column: {p_value_col}")
    
    # Filter by gene list size
    print(f"\n[2/4] Filtering by gene list size: {TEST_GENE_LIST_SIZE}")
    df_size = df[df['gene_list_size'] == TEST_GENE_LIST_SIZE].copy()
    print(f"   ✓ Filtered to {len(df_size):,} rows for size {TEST_GENE_LIST_SIZE}")
    
    # Show p-value distribution
    df_size[p_value_col] = pd.to_numeric(df_size[p_value_col], errors='coerce')
    p_values = df_size[p_value_col].dropna()
    print(f"\n   P-value statistics:")
    print(f"     Total rows: {len(p_values):,}")
    print(f"     Mean: {p_values.mean():.4f}")
    print(f"     Median: {p_values.median():.4f}")
    print(f"     Min: {p_values.min():.4f}")
    print(f"     Max: {p_values.max():.4f}")
    print(f"     Rows with p <= {TEST_P_THRESHOLD}: {(p_values <= TEST_P_THRESHOLD).sum():,} ({(p_values <= TEST_P_THRESHOLD).sum() / len(p_values) * 100:.1f}%)")
    
    # Filter by p-value threshold
    print(f"\n[3/4] Filtering by p-value threshold: <= {TEST_P_THRESHOLD}")
    df_filtered = df_size[df_size[p_value_col] <= TEST_P_THRESHOLD].copy()
    print(f"   ✓ Filtered to {len(df_filtered):,} rows with p-value <= {TEST_P_THRESHOLD}")
    
    # Filter by libraries
    if 'Library' in df_filtered.columns:
        df_filtered = df_filtered[df_filtered['Library'].isin(TEST_LIBRARIES)].copy()
        print(f"   ✓ Filtered to {len(df_filtered):,} rows for selected libraries")
    
    # Get unique permutations
    unique_permutations = df_filtered.groupby('filename').size().reset_index()
    n_permutations = len(unique_permutations)
    print(f"\n   ✓ Found {n_permutations} unique permutations")
    
    if n_permutations == 0:
        print(f"\n❌ No permutations found after filtering")
        return
    
    # Test cluster computation on a few permutations
    print(f"\n[4/4] Computing clusters for sample permutations...")
    analyzer = NetworkConnectivityAnalyzer()
    
    sample_size = min(5, n_permutations)
    sample_permutations = unique_permutations.head(sample_size)
    
    all_cluster_stats = []
    for idx, (_, perm_info) in enumerate(sample_permutations.iterrows(), 1):
        filename = perm_info['filename']
        
        # Get data for this permutation
        perm_data = df_filtered[df_filtered['filename'] == filename]
        
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
        
        # Store cluster statistics
        if clusters:
            largest_cluster = clusters[0]  # Already sorted by size
            all_cluster_stats.append({
                'filename': filename,
                'n_clusters': len(clusters),
                'largest_cluster_size': largest_cluster['size'],
                'largest_cluster_genes': largest_cluster['n_genes'],
                'largest_cluster_terms': largest_cluster['n_terms'],
                'largest_cluster_edges': largest_cluster['n_edges'],
                'largest_cluster_density': largest_cluster['density'],
            })
            print(f"   [{idx}/{sample_size}] {filename}: {len(clusters)} clusters, largest: {largest_cluster['size']} nodes")
        else:
            all_cluster_stats.append({
                'filename': filename,
                'n_clusters': 0,
                'largest_cluster_size': 0,
                'largest_cluster_genes': 0,
                'largest_cluster_terms': 0,
                'largest_cluster_edges': 0,
                'largest_cluster_density': 0.0,
            })
            print(f"   [{idx}/{sample_size}] {filename}: 0 clusters")
    
    # Show summary
    if all_cluster_stats:
        stats_df = pd.DataFrame(all_cluster_stats)
        print(f"\n" + "=" * 80)
        print("SAMPLE CLUSTER STATISTICS")
        print("=" * 80)
        print(stats_df.to_string(index=False))
        print()
        print(f"Summary for {sample_size} sample permutations:")
        print(f"  Average clusters per permutation: {stats_df['n_clusters'].mean():.2f}")
        print(f"  Average largest cluster size: {stats_df['largest_cluster_size'].mean():.2f}")
        print(f"  Average largest cluster genes: {stats_df['largest_cluster_genes'].mean():.2f}")
        print(f"  Average largest cluster terms: {stats_df['largest_cluster_terms'].mean():.2f}")
        print(f"  Average largest cluster density: {stats_df['largest_cluster_density'].mean():.4f}")
    
    # Now test the full null distribution computation
    print(f"\n" + "=" * 80)
    print("TESTING FULL NULL DISTRIBUTION COMPUTATION")
    print("=" * 80)
    
    try:
        null_stats = compute_null_distribution_from_raw_permutations(
            MERGED_PERMUTATION_FILE,
            TEST_GENE_LIST_SIZE,
            TEST_LIBRARIES,
            TEST_P_THRESHOLD
        )
        
        print(f"\n✓ Null distribution computation successful!")
        print(f"   Computed statistics for {len(null_stats)} metrics")
        print(f"\nSample statistics:")
        for metric, stats in list(null_stats.items())[:5]:
            print(f"   {metric}:")
            print(f"     Mean: {stats['mean']:.2f}")
            print(f"     Std:  {stats['std']:.2f}")
            print(f"     N:    {stats['n']}")
            print(f"     Min:  {stats['min']:.2f}")
            print(f"     Max:  {stats['max']:.2f}")
        
    except Exception as e:
        print(f"\n❌ Error computing null distribution: {e}")
        import traceback
        traceback.print_exc()
        return
    
    print(f"\n" + "=" * 80)
    print("TEST COMPLETED SUCCESSFULLY ✓")
    print("=" * 80)


if __name__ == '__main__':
    test_cluster_computation_from_permutations()
