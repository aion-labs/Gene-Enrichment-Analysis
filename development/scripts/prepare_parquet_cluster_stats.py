#!/usr/bin/env python3
"""
Prepare cluster statistics in Parquet format, segmented by gene list size.

This script:
1. Loads cluster statistics from TSV (or extracts if needed)
2. Adds boolean columns for library membership
3. Segments data by gene list size
4. Saves as Parquet files for fast loading

Output: results/permutation_cluster_statistics_parquet/cluster_stats_size_*.parquet
"""

import pandas as pd
import logging
from pathlib import Path
from typing import Set, List

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def get_all_libraries_from_data(df: pd.DataFrame) -> Set[str]:
    """Extract all unique library names from the libraries column."""
    all_libs = set()
    for libs_str in df['libraries'].dropna():
        if isinstance(libs_str, str):
            libs = [lib.strip() for lib in libs_str.split(',')]
            all_libs.update(libs)
    return all_libs


def add_library_boolean_columns(df: pd.DataFrame, all_libraries: Set[str]) -> pd.DataFrame:
    """
    Add boolean columns for each library indicating if it's present in the cluster.
    
    Args:
        df: DataFrame with 'libraries' column (comma-separated string)
        all_libraries: Set of all library names to create columns for
        
    Returns:
        DataFrame with additional boolean columns (has_<library_name>)
    """
    df = df.copy()
    
    # Create boolean columns for each library
    for lib in sorted(all_libraries):
        col_name = f"has_{lib.lower().replace(' ', '_').replace(':', '_')}"
        # Check if library name appears in the libraries string
        df[col_name] = df['libraries'].apply(
            lambda x: lib in str(x).split(',') if pd.notna(x) else False
        )
    
    return df


def prepare_parquet_files(
    input_tsv: str,
    output_dir: str,
    libraries_to_include: List[str] = None
) -> None:
    """
    Prepare Parquet files segmented by gene list size.
    
    Args:
        input_tsv: Path to cluster statistics TSV file
        output_dir: Directory to save Parquet files
        libraries_to_include: Optional list of library names to filter by
                              (if None, includes all libraries)
    """
    input_path = Path(input_tsv)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_tsv}")
    
    logger.info(f"Loading cluster statistics from {input_tsv}")
    df = pd.read_csv(input_tsv, sep='\t')
    logger.info(f"Loaded {len(df):,} cluster rows")
    
    # Get all unique libraries
    all_libraries = get_all_libraries_from_data(df)
    logger.info(f"Found {len(all_libraries)} unique libraries: {sorted(all_libraries)}")
    
    # Filter by libraries if specified
    if libraries_to_include:
        logger.info(f"Filtering to include only: {libraries_to_include}")
        # Filter clusters that contain at least one of the specified libraries
        mask = df['libraries'].apply(
            lambda x: any(lib in str(x) for lib in libraries_to_include) if pd.notna(x) else False
        )
        df = df[mask].copy()
        logger.info(f"Filtered to {len(df):,} clusters containing specified libraries")
        # Update all_libraries to only include specified ones
        all_libraries = set(libraries_to_include)
    
    # Add boolean columns for library membership
    logger.info("Adding boolean columns for library membership...")
    df = add_library_boolean_columns(df, all_libraries)
    
    # Get unique gene list sizes
    gene_list_sizes = sorted(df['gene_list_size'].unique())
    logger.info(f"Found {len(gene_list_sizes)} unique gene list sizes: {gene_list_sizes}")
    
    # Process each gene list size
    total_clusters = 0
    for size in gene_list_sizes:
        size_df = df[df['gene_list_size'] == size].copy()
        
        # Save as Parquet
        output_file = output_path / f"cluster_stats_size_{size}.parquet"
        size_df.to_parquet(output_file, index=False, compression='snappy')
        
        file_size_mb = output_file.stat().st_size / (1024 * 1024)
        logger.info(
            f"Size {size:4d}: {len(size_df):,} clusters -> {output_file.name} "
            f"({file_size_mb:.2f} MB)"
        )
        total_clusters += len(size_df)
    
    logger.info("=" * 80)
    logger.info(f"Summary:")
    logger.info(f"  Total clusters processed: {total_clusters:,}")
    logger.info(f"  Gene list sizes: {len(gene_list_sizes)}")
    logger.info(f"  Libraries: {len(all_libraries)}")
    logger.info(f"  Output directory: {output_path}")
    logger.info("=" * 80)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Prepare cluster statistics in Parquet format by gene list size'
    )
    parser.add_argument(
        '--input',
        type=str,
        default='results/permutation_cluster_statistics.tsv',
        help='Input cluster statistics TSV file'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='results/permutation_cluster_statistics_parquet',
        help='Output directory for Parquet files'
    )
    parser.add_argument(
        '--libraries',
        type=str,
        nargs='+',
        default=None,
        help='Optional: Filter to specific libraries (e.g., Reactome "GO BP")'
    )
    
    args = parser.parse_args()
    
    prepare_parquet_files(
        args.input,
        args.output_dir,
        libraries_to_include=args.libraries
    )
    
    logger.info("✓ Parquet files prepared successfully!")


if __name__ == '__main__':
    main()
