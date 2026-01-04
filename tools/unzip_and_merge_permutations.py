#!/usr/bin/env python3
"""
Unzip and merge all permutation results from zip files in the permutations folder.
Extracts gene list size from folder name and filename from file name.
Keeps only columns present in all files.
"""

import zipfile
import pandas as pd
import os
import re
from pathlib import Path
from collections import defaultdict

def extract_gene_list_size(folder_path):
    """Extract gene list size from folder name like 'size_100' -> 100"""
    match = re.search(r'size_(\d+)', folder_path)
    if match:
        return int(match.group(1))
    return None

def get_filename_without_extension(file_path):
    """Get filename without extension, e.g., 'permutation_0001.tsv' -> 'permutation_0001'"""
    return Path(file_path).stem

def find_common_columns(zip_files, sample_size=50):
    """Find columns that are present in all zip files."""
    print(f"Finding common columns by sampling {sample_size} files from each zip...")
    
    all_column_sets = []
    
    for zip_path in zip_files:
        column_sets = []
        sampled = 0
        
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            tsv_files = [f for f in zip_ref.namelist() if f.endswith('.tsv')]
            
            for tsv_file in tsv_files[:sample_size]:
                try:
                    with zip_ref.open(tsv_file) as f:
                        df = pd.read_csv(f, sep='\t', nrows=1)  # Just read header
                        column_sets.append(set(df.columns))
                        sampled += 1
                except Exception as e:
                    continue
        
        if column_sets:
            # Find intersection within this zip
            zip_common = set.intersection(*column_sets) if len(column_sets) > 1 else column_sets[0]
            all_column_sets.append(zip_common)
            print(f"  {zip_path.name}: {len(zip_common)} columns (sampled {sampled} files)")
    
    if not all_column_sets:
        raise ValueError("No files could be read to determine common columns")
    
    # Find intersection across all zips
    common_columns = set.intersection(*all_column_sets)
    print(f"\nFound {len(common_columns)} common columns across all zip files")
    print(f"Columns: {sorted(common_columns)}")
    
    return common_columns

def merge_permutation_results(zip_files, output_path, common_columns=None):
    """
    Merge all TSV files from all zip files into a single file.
    
    Args:
        zip_files: List of paths to zip files
        output_path: Path to the output merged TSV file
        common_columns: Optional set of columns to keep (if None, will be determined)
    """
    all_dataframes = []
    
    # Find common columns if not provided
    if common_columns is None:
        common_columns = find_common_columns(zip_files)
    
    total_files = 0
    
    for zip_path in zip_files:
        print(f"\nProcessing {zip_path.name}...")
        
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Get all TSV files in the zip
            tsv_files = [f for f in zip_ref.namelist() if f.endswith('.tsv')]
            
            print(f"  Found {len(tsv_files)} TSV files...")
            
            for tsv_file in tsv_files:
                # Extract folder path and filename
                folder_path = os.path.dirname(tsv_file)
                filename = os.path.basename(tsv_file)
                
                # Extract gene list size from folder name
                gene_list_size = extract_gene_list_size(folder_path)
                if gene_list_size is None:
                    continue
                
                # Get filename without extension
                filename_no_ext = get_filename_without_extension(filename)
                
                try:
                    # Read the TSV file from zip
                    with zip_ref.open(tsv_file) as f:
                        df = pd.read_csv(f, sep='\t')
                    
                    # Keep only common columns
                    available_columns = [col for col in common_columns if col in df.columns]
                    df = df[available_columns].copy()
                    
                    # Add the two new columns
                    df['gene_list_size'] = gene_list_size
                    df['filename'] = filename_no_ext
                    
                    # Reorder columns to put the new ones at the end
                    cols = [col for col in df.columns if col not in ['gene_list_size', 'filename']]
                    df = df[cols + ['gene_list_size', 'filename']]
                    
                    all_dataframes.append(df)
                    total_files += 1
                    
                    if total_files % 1000 == 0:
                        print(f"  Processed {total_files} files...")
                    
                except Exception as e:
                    print(f"  Warning: Error processing {tsv_file}: {e}")
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
    if 'Library' in merged_df.columns:
        print(f"Unique libraries: {merged_df['Library'].nunique()}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Unzip and merge all permutation results from zip files'
    )
    parser.add_argument(
        '--zip-dir',
        type=str,
        default='.',
        help='Directory containing zip files (default: current directory)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Path to output merged TSV file'
    )
    
    args = parser.parse_args()
    
    # Set default output if not provided
    if args.output is None:
        PROJECT_ROOT = Path(__file__).resolve().parent.parent
        args.output = str(PROJECT_ROOT / "sandbox" / "merged_permutation_results.tsv")
    
    # Find all zip files in the directory
    zip_dir = Path(args.zip_dir)
    zip_files = sorted(list(zip_dir.glob("*.zip")))
    
    if not zip_files:
        print(f"Error: No zip files found in {zip_dir}")
        exit(1)
    
    print(f"Found {len(zip_files)} zip files:")
    for zf in zip_files:
        print(f"  - {zf.name}")
    print()
    
    output_path = zip_dir / args.output
    
    print(f"Output will be saved to: {output_path}")
    print()
    
    # Merge the results
    merge_permutation_results(zip_files, output_path)

