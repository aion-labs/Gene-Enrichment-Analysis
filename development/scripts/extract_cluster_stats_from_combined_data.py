#!/usr/bin/env python3
"""
Extract cluster-level statistics from combined permutation data.

This script processes the combined permutation data (organized by size and library)
and extracts cluster-level statistics for each permutation's combined network.

Output columns:
- permutation_id: Unique permutation identifier (size_Permutation)
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
        permutation_id: Unique identifier for this permutation
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
    combined_data_dir: str,
    output_file: str,
    libraries_to_include: List[str] = None
) -> pd.DataFrame:
    """
    Extract cluster statistics from all permutations in combined data.
    
    Args:
        combined_data_dir: Directory containing combined permutation data (organized by size)
        output_file: Path to save cluster statistics TSV
        libraries_to_include: Optional list of library names to filter by
        
    Returns:
        DataFrame with all cluster statistics
    """
    combined_dir = Path(combined_data_dir)
    if not combined_dir.exists():
        raise FileNotFoundError(f"Combined data directory not found: {combined_data_dir}")
    
    logger.info(f"Processing combined permutation data from {combined_data_dir}")
    
    # Get all size directories
    size_dirs = sorted([d for d in combined_dir.iterdir() if d.is_dir() and d.name.startswith('size_')])
    logger.info(f"Found {len(size_dirs)} size directories")
    
    # Collect all permutation data by (size, permutation)
    # We need to combine data from all libraries for each permutation
    permutation_data = defaultdict(lambda: defaultdict(list))
    
    for size_dir in size_dirs:
        # Extract size from directory name
        try:
            gene_list_size = int(size_dir.name.replace('size_', ''))
        except ValueError:
            logger.warning(f"Skipping invalid size directory: {size_dir.name}")
            continue
        
        logger.info(f"Processing size {gene_list_size}...")
        
        # Get all library files in this size directory
        library_files = list(size_dir.glob('*.tsv'))
        logger.info(f"  Found {len(library_files)} library files")
        
        for lib_file in library_files:
            try:
                # Read library data
                lib_df = pd.read_csv(lib_file, sep='\t')
                
                # Filter by libraries if specified
                if libraries_to_include:
                    lib_df = lib_df[lib_df['Library'].isin(libraries_to_include)]
                    if lib_df.empty:
                        continue
                
                # Group by Permutation
                for perm_num, perm_group in lib_df.groupby('Permutation'):
                    permutation_id = f"{gene_list_size}_{perm_num}"
                    permutation_data[gene_list_size][perm_num].append(perm_group)
                
            except Exception as e:
                logger.warning(f"Error processing {lib_file}: {e}")
                continue
    
    logger.info(f"\nFound permutations across {len(permutation_data)} sizes")
    
    # Count total permutations
    total_perms = sum(len(perms) for perms in permutation_data.values())
    logger.info(f"Total unique permutations: {total_perms}")
    
    # Initialize analyzer (will be reused for each permutation)
    analyzer = NetworkConnectivityAnalyzer()
    
    # Process each permutation
    all_cluster_rows = []
    processed = 0
    
    for gene_list_size in sorted(permutation_data.keys()):
        for perm_num in sorted(permutation_data[gene_list_size].keys()):
            # Combine data from all libraries for this permutation
            perm_data_list = permutation_data[gene_list_size][perm_num]
            perm_data = pd.concat(perm_data_list, ignore_index=True)
            
            # Process permutation
            try:
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
                logger.warning(f"Error processing permutation {gene_list_size}_{perm_num}: {e}")
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
        description='Extract cluster statistics from combined permutation data'
    )
    parser.add_argument(
        '--combined-data-dir',
        type=str,
        default='results/combined_permutation_data',
        help='Directory containing combined permutation data (organized by size)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='results/permutation_cluster_statistics.tsv',
        help='Path to save cluster statistics TSV file'
    )
    parser.add_argument(
        '--libraries',
        type=str,
        nargs='+',
        default=None,
        help='Optional: Filter to specific libraries'
    )
    
    args = parser.parse_args()
    
    # Extract cluster statistics
    cluster_df = extract_all_cluster_statistics(
        args.combined_data_dir,
        args.output,
        libraries_to_include=args.libraries
    )
    
    logger.info(f"\n✓ Complete! Results saved to {args.output}")


if __name__ == '__main__':
    main()

