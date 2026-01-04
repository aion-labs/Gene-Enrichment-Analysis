#!/usr/bin/env python3
"""
Merge all permutation files from all folders, keeping only columns present in all files.

This script:
1. Scans all permutation result folders
2. Finds common columns across all files
3. Merges all permutation files with only common columns
4. Adds gene_list_size and filename columns
5. Saves to a single merged file
"""

import pandas as pd
import os
import re
from pathlib import Path
from collections import defaultdict
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def extract_gene_list_size(folder_path):
    """Extract gene list size from folder name like 'size_100' -> 100"""
    match = re.search(r'size_(\d+)', folder_path)
    if match:
        return int(match.group(1))
    return None


def find_all_permutation_files(base_dirs):
    """
    Find all permutation files across all base directories.
    
    Returns:
        List of tuples: (size_dir, perm_file, gene_list_size)
    """
    all_files = []
    
    for base_dir in base_dirs:
        base_path = Path(base_dir)
        if not base_path.exists():
            logger.warning(f"Directory not found: {base_dir}")
            continue
        
        # Find all size_* directories
        size_dirs = sorted([d for d in base_path.iterdir() if d.is_dir() and d.name.startswith('size_')])
        
        for size_dir in size_dirs:
            gene_list_size = extract_gene_list_size(size_dir.name)
            if gene_list_size is None:
                continue
            
            # Find all permutation files
            perm_files = sorted([f for f in size_dir.iterdir() if f.suffix == '.tsv' and f.name.startswith('permutation_')])
            
            for perm_file in perm_files:
                all_files.append((size_dir, perm_file, gene_list_size))
    
    return all_files


def find_common_columns(all_files, sample_size=100):
    """
    Find columns that are present in all permutation files.
    
    Args:
        all_files: List of (size_dir, perm_file, gene_list_size) tuples
        sample_size: Number of files to sample for finding common columns
    
    Returns:
        Set of common column names
    """
    logger.info(f"Finding common columns by sampling {sample_size} files...")
    
    # Sample files from different sizes and folders
    column_sets = []
    sampled = 0
    
    for size_dir, perm_file, gene_list_size in all_files:
        if sampled >= sample_size:
            break
        
        try:
            df = pd.read_csv(perm_file, sep='\t', nrows=1)  # Just read header
            column_sets.append(set(df.columns))
            sampled += 1
        except Exception as e:
            logger.debug(f"Error reading {perm_file}: {e}")
            continue
    
    if not column_sets:
        raise ValueError("No files could be read to determine common columns")
    
    # Find intersection of all column sets
    common_columns = set.intersection(*column_sets)
    
    logger.info(f"Sampled {sampled} files")
    logger.info(f"Found {len(common_columns)} common columns: {sorted(common_columns)}")
    
    return common_columns


def merge_all_permutation_files(
    base_dirs,
    output_file,
    common_columns=None
):
    """
    Merge all permutation files from all base directories.
    
    Args:
        base_dirs: List of base directories to search for permutation files
        output_file: Path to output merged TSV file
        common_columns: Optional set of columns to keep (if None, will be determined)
    """
    logger.info("="*80)
    logger.info("Merging All Permutation Files")
    logger.info("="*80)
    
    # Find all permutation files
    logger.info("Scanning for permutation files...")
    all_files = find_all_permutation_files(base_dirs)
    logger.info(f"Found {len(all_files)} permutation files to process")
    
    if not all_files:
        raise ValueError("No permutation files found!")
    
    # Find common columns if not provided
    if common_columns is None:
        common_columns = find_common_columns(all_files)
    
    # Always add these columns (will be added during processing)
    columns_to_keep = list(common_columns)
    logger.info(f"Keeping {len(columns_to_keep)} common columns")
    
    # Process all files
    all_dataframes = []
    processed = 0
    errors = 0
    
    for size_dir, perm_file, gene_list_size in all_files:
        try:
            # Read file
            df = pd.read_csv(perm_file, sep='\t')
            
            # Keep only common columns
            available_columns = [col for col in columns_to_keep if col in df.columns]
            df = df[available_columns].copy()
            
            # Add gene_list_size and filename
            df['gene_list_size'] = gene_list_size
            df['filename'] = perm_file.stem  # e.g., "permutation_0001"
            
            # Reorder columns: put new columns at the end
            cols = [col for col in df.columns if col not in ['gene_list_size', 'filename']]
            df = df[cols + ['gene_list_size', 'filename']]
            
            all_dataframes.append(df)
            processed += 1
            
            # Progress logging
            if processed % 1000 == 0:
                logger.info(f"Processed {processed}/{len(all_files)} files...")
                
        except Exception as e:
            logger.warning(f"Error processing {perm_file}: {e}")
            errors += 1
            continue
    
    logger.info(f"\nProcessed {processed} files successfully, {errors} errors")
    
    if not all_dataframes:
        raise ValueError("No dataframes to merge!")
    
    # Concatenate all dataframes
    logger.info("Merging all dataframes...")
    merged_df = pd.concat(all_dataframes, ignore_index=True)
    
    # Save to output file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Saving merged results to {output_file}...")
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    # Print summary
    logger.info("\n" + "="*80)
    logger.info("SUMMARY")
    logger.info("="*80)
    logger.info(f"Total files processed: {processed}")
    logger.info(f"Total rows: {len(merged_df):,}")
    logger.info(f"Columns: {len(merged_df.columns)}")
    logger.info(f"Columns kept: {', '.join(sorted(merged_df.columns))}")
    logger.info(f"\nGene list sizes included: {sorted(merged_df['gene_list_size'].unique())}")
    logger.info(f"Unique permutations: {merged_df['filename'].nunique():,}")
    if 'Library' in merged_df.columns:
        logger.info(f"Unique libraries: {merged_df['Library'].nunique()}")
        logger.info(f"Libraries: {', '.join(sorted(merged_df['Library'].unique()))}")
    
    # Count by size
    logger.info(f"\nRows by gene list size:")
    size_counts = merged_df.groupby('gene_list_size').size()
    for size, count in size_counts.items():
        n_perms = merged_df[merged_df['gene_list_size'] == size]['filename'].nunique()
        logger.info(f"  Size {size:4d}: {count:8,} rows from {n_perms:4d} permutations")
    
    logger.info("="*80)
    
    return merged_df


def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Merge all permutation files from all folders, keeping only common columns'
    )
    parser.add_argument(
        '--base-dirs',
        type=str,
        nargs='+',
        default=[
            'results/permutations_results-50-500',
            'results/permutation_results-550-1000',
            'results/permutation_results-FirstRun-50-to-1000'
        ],
        help='Base directories containing permutation results (size_* folders)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='results/merged_permutation_results_all.tsv',
        help='Path to output merged TSV file'
    )
    
    args = parser.parse_args()
    
    logger.info(f"Base directories to process: {args.base_dirs}")
    logger.info(f"Output file: {args.output}")
    logger.info("")
    
    # Merge all files
    merged_df = merge_all_permutation_files(
        args.base_dirs,
        args.output
    )
    
    logger.info(f"\n✓ Complete! Results saved to {args.output}")


if __name__ == '__main__':
    main()

