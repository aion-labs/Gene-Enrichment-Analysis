#!/usr/bin/env python3
"""
Extract cluster-level statistics from all permutation files.

This script processes all permutations and extracts cluster-level statistics
for each cluster in each permutation's combined network. Each row represents
one cluster from one permutation.

Output columns:
- filename: Permutation filename
- gene_list_size: Size of the gene list (50, 100, 150, etc.)
- cluster_number: Cluster rank by size (1 = largest)
- cluster_size: Total nodes (genes + terms) in cluster
- n_genes: Number of genes in cluster
- n_terms: Number of terms in cluster
- n_edges: Number of edges in cluster
- n_libraries: Number of unique libraries in cluster
- density: Cluster density (edges / (genes × libraries))
- libraries: Comma-separated list of libraries in cluster
"""

import pandas as pd
import logging
import sys
from pathlib import Path
from typing import List, Dict

# Add code directory to path to import code modules
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))
from network_connectivity_benchmark import NetworkConnectivityAnalyzer

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_permutation(
    permutation_data: pd.DataFrame,
    filename: str,
    gene_list_size: int,
    analyzer: NetworkConnectivityAnalyzer
) -> List[Dict]:
    """
    Process a single permutation and extract all cluster statistics.
    
    Args:
        permutation_data: DataFrame with iGEA results for this permutation
        filename: Permutation filename
        gene_list_size: Size of the gene list
        analyzer: NetworkConnectivityAnalyzer instance (will be reset)
        
    Returns:
        List of dictionaries, one per cluster
    """
    # Reset analyzer for this permutation
    analyzer.reset()
    
    # Convert permutation data to iGEA results format
    results = []
    for _, row in permutation_data.iterrows():
        # Parse genes removed
        genes_removed = row.get('Genes removed for next iteration', '')
        if isinstance(genes_removed, str):
            genes_list = [g.strip() for g in genes_removed.split(',') if g.strip()]
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
    
    # Get all clusters
    clusters = analyzer.get_clusters()
    
    # Extract cluster statistics
    cluster_rows = []
    for cluster in clusters:
        # Get libraries in this cluster
        libraries_in_cluster = set()
        for term in cluster['terms']:
            library = analyzer.term_to_library.get(term, "Unknown")
            libraries_in_cluster.add(library)
        
        libraries_str = ', '.join(sorted(libraries_in_cluster))
        
        cluster_rows.append({
            'filename': filename,
            'gene_list_size': gene_list_size,
            'cluster_number': cluster['cluster_number'],
            'cluster_size': cluster['size'],
            'n_genes': cluster['n_genes'],
            'n_terms': cluster['n_terms'],
            'n_edges': cluster['n_edges'],
            'n_libraries': cluster['n_libraries'],
            'density': cluster['density'],
            'libraries': libraries_str,
        })
    
    return cluster_rows


def extract_all_cluster_statistics(
    merged_permutation_file: str,
    output_file: str
) -> pd.DataFrame:
    """
    Extract cluster statistics from all permutations.
    
    Args:
        merged_permutation_file: Path to merged permutation results TSV
        output_file: Path to save cluster statistics TSV
        
    Returns:
        DataFrame with all cluster statistics
    """
    logger.info(f"Loading permutation results from {merged_permutation_file}")
    df = pd.read_csv(merged_permutation_file, sep='\t')
    logger.info(f"Loaded {len(df)} permutation result rows")
    
    # Get unique permutations
    unique_permutations = df.groupby(['filename', 'gene_list_size']).size().reset_index()
    n_permutations = len(unique_permutations)
    logger.info(f"Found {n_permutations} unique permutations")
    
    # Initialize analyzer (will be reused for each permutation)
    analyzer = NetworkConnectivityAnalyzer()
    
    # Process each permutation
    all_cluster_rows = []
    for idx, (_, perm_info) in enumerate(unique_permutations.iterrows(), 1):
        filename = perm_info['filename']
        gene_list_size = perm_info['gene_list_size']
        
        # Get data for this permutation
        perm_data = df[(df['filename'] == filename) & 
                       (df['gene_list_size'] == gene_list_size)]
        
        # Process permutation
        cluster_rows = process_permutation(
            perm_data,
            filename,
            gene_list_size,
            analyzer
        )
        
        all_cluster_rows.extend(cluster_rows)
        
        # Progress logging
        if idx % 100 == 0:
            logger.info(f"Processed {idx}/{n_permutations} permutations "
                       f"({len(all_cluster_rows)} clusters so far)...")
    
    # Create DataFrame
    cluster_df = pd.DataFrame(all_cluster_rows)
    
    # Sort by gene_list_size, then filename, then cluster_number
    cluster_df = cluster_df.sort_values(
        ['gene_list_size', 'filename', 'cluster_number']
    ).reset_index(drop=True)
    
    # Save to file
    cluster_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved {len(cluster_df)} cluster statistics to {output_file}")
    
    # Print summary
    logger.info("\n" + "="*80)
    logger.info("SUMMARY")
    logger.info("="*80)
    logger.info(f"Total permutations processed: {n_permutations}")
    logger.info(f"Total clusters: {len(cluster_df)}")
    logger.info(f"Average clusters per permutation: {len(cluster_df) / n_permutations:.2f}")
    logger.info(f"\nClusters by gene list size:")
    size_counts = cluster_df.groupby('gene_list_size').size()
    for size, count in size_counts.items():
        n_perms = len(cluster_df[cluster_df['gene_list_size'] == size]['filename'].unique())
        logger.info(f"  Size {size}: {count} clusters from {n_perms} permutations "
                   f"(avg {count/n_perms:.2f} clusters/permutation)")
    
    return cluster_df


def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Extract cluster statistics from all permutation files'
    )
    parser.add_argument(
        '--permutation-file',
        type=str,
        default='archive/permutation_analysis/results/permutation_results/merged_permutation_results.tsv',
        help='Path to merged permutation results TSV file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Path to save cluster statistics TSV file'
    )
    
    args = parser.parse_args()
    
    # Set default output if not provided
    if args.output is None:
        PROJECT_ROOT = Path(__file__).resolve().parent.parent
        args.output = str(PROJECT_ROOT / "sandbox" / "permutation_cluster_statistics.tsv")
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Extract cluster statistics
    cluster_df = extract_all_cluster_statistics(
        args.permutation_file,
        args.output
    )
    
    logger.info(f"\n✓ Complete! Results saved to {args.output}")


if __name__ == '__main__':
    main()
