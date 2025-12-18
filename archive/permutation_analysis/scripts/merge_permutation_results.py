#!/usr/bin/env python3
"""
Script to merge all permutation results from a zip file.
Extracts gene list size from folder name and filename from file name.
Skips the Description column.
"""

import zipfile
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

def get_filename_without_extension(file_path):
    """Get filename without extension, e.g., 'permutation_0001.tsv' -> 'permutation_0001'"""
    return Path(file_path).stem

def merge_permutation_results(zip_path, output_path):
    """
    Merge all TSV files from the zip into a single file.
    
    Args:
        zip_path: Path to the zip file containing permutation results
        output_path: Path to the output merged TSV file
    """
    all_dataframes = []
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Get all TSV files in the zip
        tsv_files = [f for f in zip_ref.namelist() if f.endswith('.tsv')]
        
        print(f"Found {len(tsv_files)} TSV files to process...")
        
        for tsv_file in tsv_files:
            # Extract folder path and filename
            folder_path = os.path.dirname(tsv_file)
            filename = os.path.basename(tsv_file)
            
            # Extract gene list size from folder name
            gene_list_size = extract_gene_list_size(folder_path)
            if gene_list_size is None:
                print(f"Warning: Could not extract gene list size from {folder_path}, skipping...")
                continue
            
            # Get filename without extension
            filename_no_ext = get_filename_without_extension(filename)
            
            try:
                # Read the TSV file from zip
                with zip_ref.open(tsv_file) as f:
                    df = pd.read_csv(f, sep='\t')
                
                # Skip Description column if it exists
                if 'Description' in df.columns:
                    df = df.drop(columns=['Description'])
                
                # Add the two new columns
                df['gene_list_size'] = gene_list_size
                df['filename'] = filename_no_ext
                
                # Reorder columns to put the new ones at the end (or beginning if preferred)
                # Move new columns to the end
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
    print("Merging all dataframes...")
    merged_df = pd.concat(all_dataframes, ignore_index=True)
    
    # Save to output file
    print(f"Saving merged results to {output_path}...")
    merged_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"Done! Merged {len(all_dataframes)} files into {output_path}")
    print(f"Total rows: {len(merged_df)}")
    print(f"Columns: {', '.join(merged_df.columns)}")

if __name__ == "__main__":
    import sys
    
    # Set paths
    script_dir = Path(__file__).parent
    permutation_results_dir = script_dir / "permutation_results"
    
    # Find zip file(s) in the permutation_results directory
    zip_files = list(permutation_results_dir.glob("*.zip"))
    
    if not zip_files:
        print(f"Error: No zip files found in {permutation_results_dir}")
        sys.exit(1)
    
    # Use the first zip file found (or you can modify to handle multiple)
    zip_path = zip_files[0]
    if len(zip_files) > 1:
        print(f"Warning: Multiple zip files found. Using: {zip_path.name}")
    
    output_path = permutation_results_dir / "merged_permutation_results.tsv"
    
    # Check if zip file exists
    if not zip_path.exists():
        print(f"Error: Zip file not found at {zip_path}")
        sys.exit(1)
    
    print(f"Processing zip file: {zip_path}")
    print(f"Output will be saved to: {output_path}")
    print()
    
    # Merge the results
    merge_permutation_results(zip_path, output_path)

