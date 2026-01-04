#!/usr/bin/env python3
"""
Extract cluster-level statistics from FirstRun permutation data.

This script processes the FirstRun permutation files directly and extracts
cluster-level statistics for each permutation's combined network.
"""

import pandas as pd
import logging
from pathlib import Path
from typing import List, Dict
from collections import defaultdict
from code.network_connectivity_benchmark import NetworkConnectivityAnalyzer

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_permutation(
    permutation_data: pd.DataFrame,
    perm_num: int,
    gene_list_size: int,
    analyzer: NetworkConnectivityAnalyzer
) -> List[Dict]:
    """
    Process a single permutation and extract all cluster statistics.
    
    Args:
        permutation_data: DataFrame with iGEA results for this permutation
        perm_num: Permutation number
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
            'filename': f"permutation_{perm_num:04d}",
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
    firstrun_dir: str,
    output_file: str
) -> pd.DataFrame:
    """
    Extract cluster statistics from all permutations in FirstRun folder.
    
    Args:
        firstrun_dir: Directory containing FirstRun permutation data
        output_file: Path to save cluster statistics TSV
        
    Returns:
        DataFrame with all cluster statistics
    """
    firstrun_path = Path(firstrun_dir)
    if not firstrun_path.exists():
        raise FileNotFoundError(f"FirstRun directory not found: {firstrun_dir}")
    
    logger.info(f"Processing FirstRun permutation data from {firstrun_dir}")
    
    # Get all size directories
    size_dirs = sorted([d for d in firstrun_path.iterdir() if d.is_dir() and d.name.startswith('size_')])
    logger.info(f"Found {len(size_dirs)} size directories")
    
    # Collect all permutation files
    permutation_files = []
    for size_dir in size_dirs:
        try:
            gene_list_size = int(size_dir.name.replace('size_', ''))
        except ValueError:
            logger.warning(f"Skipping invalid size directory: {size_dir.name}")
            continue
        
        # Get all permutation files in this size directory
        perm_files = sorted(size_dir.glob('permutation_*.tsv'))
        for perm_file in perm_files:
            try:
                # Extract permutation number from filename
                perm_num = int(perm_file.stem.split('_')[1])
                permutation_files.append((gene_list_size, perm_num, perm_file))
            except (ValueError, IndexError):
                logger.warning(f"Could not parse permutation number from {perm_file}")
                continue
    
    total_perms = len(permutation_files)
    logger.info(f"Found {total_perms} permutation files to process")
    
    # Initialize analyzer (will be reused for each permutation)
    analyzer = NetworkConnectivityAnalyzer()
    
    # Process each permutation
    all_cluster_rows = []
    processed = 0
    
    for gene_list_size, perm_num, perm_file in permutation_files:
        try:
            # Read permutation file
            perm_data = pd.read_csv(perm_file, sep='\t')
            
            # Process permutation
            cluster_rows = process_permutation(
                perm_data,
                perm_num,
                gene_list_size,
                analyzer
            )
            all_cluster_rows.extend(cluster_rows)
            processed += 1
            
            # Progress logging
            if processed % 100 == 0:
                logger.info(f"Processed {processed}/{total_perms} permutations "
                           f"({len(all_cluster_rows)} clusters so far)...")
        except Exception as e:
            logger.warning(f"Error processing {perm_file}: {e}")
            continue
    
    # Create DataFrame
    cluster_df = pd.DataFrame(all_cluster_rows)
    
    if cluster_df.empty:
        logger.warning("No cluster statistics extracted!")
        return cluster_df
    
    # Sort by gene_list_size, then filename, then cluster_number
    cluster_df = cluster_df.sort_values(
        ['gene_list_size', 'filename', 'cluster_number']
    ).reset_index(drop=True)
    
    # Save to file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved {len(cluster_df)} cluster statistics to {output_file}")
    
    # Print summary
    logger.info("\n" + "="*80)
    logger.info("SUMMARY")
    logger.info("="*80)
    logger.info(f"Total permutations processed: {processed}")
    logger.info(f"Total clusters: {len(cluster_df)}")
    if processed > 0:
        logger.info(f"Average clusters per permutation: {len(cluster_df) / processed:.2f}")
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
        description='Extract cluster statistics from FirstRun permutation data'
    )
    parser.add_argument(
        '--firstrun-dir',
        type=str,
        default='results/permutation_results-FirstRun-50-to-1000',
        help='Directory containing FirstRun permutation data'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='results/permutation_cluster_statistics.tsv',
        help='Path to save cluster statistics TSV file'
    )
    
    args = parser.parse_args()
    
    # Extract cluster statistics
    cluster_df = extract_all_cluster_statistics(
        args.firstrun_dir,
        args.output
    )
    
    logger.info(f"\n✓ Complete! Results saved to {args.output}")


if __name__ == '__main__':
    main()

