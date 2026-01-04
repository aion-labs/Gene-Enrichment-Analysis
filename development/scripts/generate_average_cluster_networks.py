#!/usr/bin/env python3
"""
Find clusters matching average clustering parameters for size 300,
randomly select 3, and generate DOT graph networks.

Terms are color coded by libraries, genes are white.
"""

import pandas as pd
import numpy as np
import random
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import sys

# Add code directory to path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from network_connectivity_benchmark import NetworkConnectivityAnalyzer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_cluster_stats_for_size(parquet_dir: Path, gene_list_size: int, largest_only: bool = True) -> pd.DataFrame:
    """
    Load cluster statistics for a specific gene list size.
    
    Args:
        parquet_dir: Directory containing parquet files
        gene_list_size: Gene list size to load
        largest_only: If True, only load largest cluster (cluster_number == 1) from each permutation
    
    Returns:
        DataFrame with cluster statistics
    """
    parquet_file = parquet_dir / f"cluster_stats_size_{gene_list_size}.parquet"
    if not parquet_file.exists():
        raise FileNotFoundError(f"Parquet file not found: {parquet_file}")
    
    df = pd.read_parquet(parquet_file)
    logger.info(f"Loaded {len(df):,} clusters from size {gene_list_size}")
    
    if largest_only:
        # Filter to only largest cluster (cluster_number == 1) from each permutation
        df = df[df['cluster_number'] == 1].copy()
        logger.info(f"Filtered to {len(df):,} largest clusters (one per permutation)")
    
    return df


def calculate_average_parameters(df: pd.DataFrame) -> Dict[str, float]:
    """
    Calculate average clustering parameters from the dataframe.
    
    Since we're only looking at largest clusters (one per permutation),
    we can directly calculate means without grouping by filename.
    """
    # Since df already contains only largest clusters (one per permutation),
    # we can directly calculate means
    avg_params = {
        'avg_cluster_size': df['cluster_size'].mean(),
        'avg_n_genes': df['n_genes'].mean(),
        'avg_n_terms': df['n_terms'].mean(),
        'avg_n_edges': df['n_edges'].mean(),
        'avg_density': df['density'].mean(),
        'avg_n_libraries': df['n_libraries'].mean(),
    }
    
    logger.info("Average clustering parameters (largest clusters only):")
    for key, value in avg_params.items():
        logger.info(f"  {key}: {value:.4f}")
    
    logger.info(f"\nStatistics based on {len(df)} largest clusters")
    
    return avg_params


def find_matching_clusters(df: pd.DataFrame, avg_params: Dict[str, float], 
                          threshold: float = 0.2) -> pd.DataFrame:
    """
    Find clusters that match closely to average parameters.
    
    Args:
        df: DataFrame with cluster statistics
        avg_params: Dictionary of average parameter values
        threshold: Fractional deviation allowed (e.g., 0.2 = 20%)
    
    Returns:
        DataFrame with matching clusters
    """
    matching_clusters = []
    
    for idx, row in df.iterrows():
        # Calculate deviations for each parameter
        deviations = []
        
        # Check cluster size
        if avg_params['avg_cluster_size'] > 0:
            size_dev = abs(row['cluster_size'] - avg_params['avg_cluster_size']) / avg_params['avg_cluster_size']
            deviations.append(size_dev)
        
        # Check n_genes
        if avg_params['avg_n_genes'] > 0:
            genes_dev = abs(row['n_genes'] - avg_params['avg_n_genes']) / avg_params['avg_n_genes']
            deviations.append(genes_dev)
        
        # Check n_terms
        if avg_params['avg_n_terms'] > 0:
            terms_dev = abs(row['n_terms'] - avg_params['avg_n_terms']) / avg_params['avg_n_terms']
            deviations.append(terms_dev)
        
        # Check density (avg_edges_per_gene)
        if avg_params['avg_density'] > 0:
            density_dev = abs(row['density'] - avg_params['avg_density']) / avg_params['avg_density']
            deviations.append(density_dev)
        
        # If all deviations are within threshold, include this cluster
        if deviations and max(deviations) <= threshold:
            matching_clusters.append(idx)
    
    matching_df = df.loc[matching_clusters].copy()
    logger.info(f"Found {len(matching_df)} clusters matching average parameters (threshold: {threshold*100}%)")
    
    return matching_df


def load_permutation_data(merged_permutation_file: Path, filename: str, 
                         gene_list_size: int) -> pd.DataFrame:
    """Load permutation data for a specific filename and gene list size."""
    logger.info(f"Loading permutation data for {filename}, size {gene_list_size}")
    
    # Read the merged permutation file
    df = pd.read_csv(merged_permutation_file, sep='\t')
    
    # Filter by filename and gene list size
    perm_data = df[(df['filename'] == filename) & 
                   (df['gene_list_size'] == gene_list_size)].copy()
    
    logger.info(f"Loaded {len(perm_data)} rows for {filename}")
    return perm_data


def reconstruct_network_from_permutation(perm_data: pd.DataFrame) -> NetworkConnectivityAnalyzer:
    """Reconstruct network from permutation data."""
    analyzer = NetworkConnectivityAnalyzer()
    
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
        
        term = row.get('Term', '')
        if term:
            results.append({
                'Term': term,
                'Iteration': row.get('Iteration', 1),
                'Library': row.get('Library', ''),
                'Genes removed for next iteration': genes_list,
            })
    
    # Add to analyzer
    analyzer.add_igea_results(results)
    
    return analyzer


def generate_colored_dot_graph(analyzer: NetworkConnectivityAnalyzer, 
                               cluster: Dict,
                               output_file: Path,
                               library_colors: Dict[str, str] = None):
    """
    Generate DOT graph with color coding:
    - Terms: colored by library
    - Genes: white
    
    Args:
        analyzer: NetworkConnectivityAnalyzer with the network
        cluster: Cluster dictionary with terms list
        output_file: Path to save DOT file
        library_colors: Optional dictionary mapping library names to colors
    """
    # Default color palette for libraries
    if library_colors is None:
        # Generate colors for libraries
        libraries = set()
        for term in cluster['terms']:
            library = analyzer.term_to_library.get(term, "Unknown")
            libraries.add(library)
        
        # Color palette (distinct colors)
        color_palette = [
            '#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8',
            '#F7DC6F', '#BB8FCE', '#85C1E2', '#F8B739', '#52BE80',
            '#EC7063', '#5DADE2', '#F39C12', '#16A085', '#8E44AD'
        ]
        
        library_colors = {}
        for idx, lib in enumerate(sorted(libraries)):
            library_colors[lib] = color_palette[idx % len(color_palette)]
    
    # Get terms and genes in this cluster
    cluster_terms = set(cluster['terms'])
    cluster_genes = set()
    
    for term in cluster_terms:
        genes = analyzer.term_to_genes.get(term, set())
        cluster_genes.update(genes)
    
    # Generate DOT content
    dot_lines = ['graph cluster_network {']
    dot_lines.append('  graph [layout=neato, overlap=false, splines=true];')
    dot_lines.append('  node [shape=ellipse, style=filled];')
    dot_lines.append('')
    
    # Add term nodes (colored by library)
    for term in sorted(cluster_terms):
        library = analyzer.term_to_library.get(term, "Unknown")
        color = library_colors.get(library, '#CCCCCC')
        # Escape special characters in term name
        term_escaped = term.replace('"', '\\"').replace('\n', ' ')
        dot_lines.append(f'  "{term_escaped}" [fillcolor="{color}", label="{term_escaped}"];')
    
    # Add gene nodes (white)
    for gene in sorted(cluster_genes):
        gene_escaped = gene.replace('"', '\\"').replace('\n', ' ')
        dot_lines.append(f'  "{gene_escaped}" [fillcolor="white", label="{gene_escaped}"];')
    
    dot_lines.append('')
    
    # Add edges (gene-term connections)
    edge_count = 0
    for term in cluster_terms:
        genes = analyzer.term_to_genes.get(term, set())
        term_escaped = term.replace('"', '\\"').replace('\n', ' ')
        for gene in genes:
            if gene in cluster_genes:  # Only add edges for genes in this cluster
                gene_escaped = gene.replace('"', '\\"').replace('\n', ' ')
                dot_lines.append(f'  "{gene_escaped}" -- "{term_escaped}";')
                edge_count += 1
    
    dot_lines.append('}')
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(dot_lines))
    
    logger.info(f"Saved DOT graph to {output_file}")
    logger.info(f"  Terms: {len(cluster_terms)}, Genes: {len(cluster_genes)}, Edges: {edge_count}")


def main():
    """Main function."""
    # Configuration
    GENE_LIST_SIZE = 300
    THRESHOLD = 0.2  # 20% deviation allowed
    N_SELECT = 3  # Number of clusters to randomly select
    
    # Paths
    PROJECT_ROOT = Path(__file__).resolve().parent
    parquet_dir = PROJECT_ROOT / "permutations" / "permutation_cluster_statistics_parquet"
    merged_permutation_file = PROJECT_ROOT / "permutations" / "merged_permutation_results.tsv"
    output_dir = PROJECT_ROOT / "results" / "average_cluster_networks"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("Finding Average Clusters and Generating Network Graphs")
    logger.info("=" * 80)
    logger.info(f"Gene list size: {GENE_LIST_SIZE}")
    logger.info(f"Matching threshold: {THRESHOLD*100}%")
    logger.info(f"Number to select: {N_SELECT}")
    logger.info("")
    
    # Step 1: Load cluster statistics for size 300 (largest clusters only)
    logger.info("Step 1: Loading cluster statistics (largest clusters only)...")
    df = load_cluster_stats_for_size(parquet_dir, GENE_LIST_SIZE, largest_only=True)
    
    # Step 2: Calculate average parameters (from largest clusters only)
    logger.info("\nStep 2: Calculating average clustering parameters (largest clusters only)...")
    avg_params = calculate_average_parameters(df)
    
    # Step 3: Find matching largest clusters
    logger.info(f"\nStep 3: Finding largest clusters matching average parameters (threshold: {THRESHOLD*100}%)...")
    matching_df = find_matching_clusters(df, avg_params, threshold=THRESHOLD)
    
    if len(matching_df) == 0:
        logger.error("No matching clusters found! Try increasing the threshold.")
        return
    
    # Step 4: Randomly select N clusters
    logger.info(f"\nStep 4: Randomly selecting {N_SELECT} clusters...")
    if len(matching_df) < N_SELECT:
        logger.warning(f"Only {len(matching_df)} matching clusters found, selecting all of them")
        selected_clusters = matching_df
    else:
        selected_clusters = matching_df.sample(n=N_SELECT, random_state=42)
    
    logger.info(f"Selected {len(selected_clusters)} clusters:")
    for idx, (_, row) in enumerate(selected_clusters.iterrows(), 1):
        logger.info(f"  {idx}. {row['filename']} - Cluster {row['cluster_number']}: "
                   f"{row['cluster_size']} nodes, {row['n_genes']} genes, {row['n_terms']} terms")
    
    # Step 5: Generate DOT graphs for each selected cluster
    logger.info(f"\nStep 5: Generating DOT graphs...")
    
    # Load merged permutation file once
    logger.info(f"Loading merged permutation file: {merged_permutation_file}")
    all_perm_data = pd.read_csv(merged_permutation_file, sep='\t')
    logger.info(f"Loaded {len(all_perm_data):,} permutation result rows")
    
    for idx, (_, cluster_row) in enumerate(selected_clusters.iterrows(), 1):
        filename = cluster_row['filename']
        cluster_num = cluster_row['cluster_number']
        
        logger.info(f"\nProcessing cluster {idx}/{len(selected_clusters)}: {filename}, cluster {cluster_num}")
        
        # Load permutation data for this filename
        perm_data = all_perm_data[
            (all_perm_data['filename'] == filename) & 
            (all_perm_data['gene_list_size'] == GENE_LIST_SIZE)
        ].copy()
        
        if len(perm_data) == 0:
            logger.warning(f"No permutation data found for {filename}, skipping")
            continue
        
        # Reconstruct network
        analyzer = reconstruct_network_from_permutation(perm_data)
        
        # Get all clusters to find the specific one
        all_clusters = analyzer.get_clusters()
        
        # Find the cluster matching cluster_number
        target_cluster = None
        if cluster_num <= len(all_clusters):
            # Clusters are sorted by size, but we need to match by cluster_number
            # The cluster_number in parquet might not match the index
            # Let's try to match by size and other metrics
            for cluster in all_clusters:
                if (abs(cluster['size'] - cluster_row['cluster_size']) < 2 and
                    abs(cluster['n_genes'] - cluster_row['n_genes']) < 2 and
                    abs(cluster['n_terms'] - cluster_row['n_terms']) < 2):
                    target_cluster = cluster
                    break
            
            # Fallback: use cluster at index (cluster_num - 1) if sorted by size
            if target_cluster is None and cluster_num <= len(all_clusters):
                # Try to get by index (clusters are sorted by size, largest first)
                target_cluster = all_clusters[cluster_num - 1]
        
        if target_cluster is None:
            logger.warning(f"Could not find matching cluster for {filename}, cluster {cluster_num}")
            # Use largest cluster as fallback
            if all_clusters:
                target_cluster = all_clusters[0]
                logger.info(f"Using largest cluster as fallback: {target_cluster['size']} nodes")
            else:
                logger.error(f"No clusters found for {filename}, skipping")
                continue
        
        # Generate DOT graph
        output_file = output_dir / f"average_cluster_{idx}_{filename}_cluster{cluster_num}.dot"
        generate_colored_dot_graph(analyzer, target_cluster, output_file)
    
    logger.info("\n" + "=" * 80)
    logger.info("Done! Generated DOT graphs:")
    for idx in range(1, len(selected_clusters) + 1):
        logger.info(f"  {output_dir / f'average_cluster_{idx}_*.dot'}")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()

