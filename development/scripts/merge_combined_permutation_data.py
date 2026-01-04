#!/usr/bin/env python3
"""
Script to merge all permutation results from combined_permutation_data folder.
Extracts gene list size from folder name and adds it as a column.
Works with the combined_permutation_data structure (one file per library per size).
"""

import pandas as pd
import os
import re
from pathlib import Path

def extract_gene_list_size(folder_path):
    """Extract gene list size from folder name like 'size_100' -> 100"""
    match = re.search(r'size_(\d+)', folder_path)
    if match:
        return int(match.group(1))
    return None

def merge_combined_permutation_data(input_dir, output_path):
    """
    Merge all TSV files from the combined_permutation_data folder into a single file.
    
    Args:
        input_dir: Path to the combined_permutation_data directory
        output_path: Path to the output merged TSV file
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    all_dataframes = []
    
    # Get all size directories
    size_dirs = sorted([d for d in input_path.iterdir() if d.is_dir() and d.name.startswith('size_')])
    
    print(f"Found {len(size_dirs)} size directories to process...")
    
    for size_dir in size_dirs:
        # Extract gene list size from folder name
        gene_list_size = extract_gene_list_size(size_dir.name)
        if gene_list_size is None:
            print(f"Warning: Could not extract gene list size from {size_dir.name}, skipping...")
            continue
        
        print(f"Processing size {gene_list_size}...")
        
        # Get all TSV files in this size directory
        tsv_files = sorted([f for f in size_dir.iterdir() if f.suffix == '.tsv'])
        
        for tsv_file in tsv_files:
            try:
                # Read the TSV file
                df = pd.read_csv(tsv_file, sep='\t')
                
                # Skip Description column if it exists
                if 'Description' in df.columns:
                    df = df.drop(columns=['Description'])
                
                # Add gene_list_size column
                df['gene_list_size'] = gene_list_size
                
                # Create filename from Permutation column if it exists
                # Format: permutation_XXXX where XXXX is the permutation number
                if 'Permutation' in df.columns:
                    df['filename'] = df['Permutation'].apply(
                        lambda x: f"permutation_{int(x):04d}" if pd.notna(x) else "unknown"
                    )
                else:
                    # If no Permutation column, use a default
                    df['filename'] = f"permutation_unknown"
                
                # Reorder columns to put new ones at the end
                cols = [col for col in df.columns if col not in ['gene_list_size', 'filename']]
                df = df[cols + ['gene_list_size', 'filename']]
                
                all_dataframes.append(df)
                
            except Exception as e:
                print(f"Error processing {tsv_file}: {e}")
                continue
    
    if not all_dataframes:
        print("No dataframes to merge!")
        return
    
    # Concatenate all dataframes
    print(f"\nMerging {len(all_dataframes)} files...")
    merged_df = pd.concat(all_dataframes, ignore_index=True)
    
    # Save to output file
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Saving merged results to {output_path}...")
    merged_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"\nDone! Merged {len(all_dataframes)} files into {output_path}")
    print(f"Total rows: {len(merged_df):,}")
    print(f"Columns: {', '.join(merged_df.columns)}")
    print(f"\nGene list sizes included: {sorted(merged_df['gene_list_size'].unique())}")
    print(f"Unique permutations: {merged_df['filename'].nunique():,}")
    print(f"Unique libraries: {merged_df['Library'].nunique() if 'Library' in merged_df.columns else 'N/A'}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Merge all permutation results from combined_permutation_data folder'
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        default='results/combined_permutation_data',
        help='Path to combined_permutation_data directory'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='results/merged_permutation_results.tsv',
        help='Path to output merged TSV file'
    )
    
    args = parser.parse_args()
    
    print(f"Processing directory: {args.input_dir}")
    print(f"Output will be saved to: {args.output}")
    print()
    
    # Merge the results
    merge_combined_permutation_data(args.input_dir, args.output)

