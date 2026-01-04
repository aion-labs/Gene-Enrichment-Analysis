#!/usr/bin/env python3
"""
Combine permutation data from multiple runs into a single organized structure.
Organizes by gene list size, with one combined file per library.
"""

import pandas as pd
import os
from pathlib import Path
from collections import defaultdict
import json

def get_libraries_in_size(size_dir):
    """Get unique libraries from a size directory by checking multiple permutation files.
    
    Checks multiple files because a library might not have significant results
    in all permutations, but we want to know if it exists at this size at all.
    """
    if not os.path.exists(size_dir):
        return set()
    
    libraries = set()
    # Check multiple permutation files to get all libraries that exist at this size
    # (some libraries might not have results in all permutations)
    perm_files = sorted([f for f in os.listdir(size_dir) if f.endswith('.tsv')])
    
    # Check first 10 files (or all if fewer than 10)
    files_to_check = perm_files[:min(10, len(perm_files))]
    
    for perm_file in files_to_check:
        file_path = os.path.join(size_dir, perm_file)
        try:
            df = pd.read_csv(file_path, sep='\t')
            libraries.update(df['Library'].unique())
        except Exception as e:
            continue
    
    return libraries

def find_libraries_with_high_coverage(min_coverage=0.85):
    """Find libraries with high coverage (>= min_coverage) across sizes from 50 to 1000.
    
    Checks all three permutation result folders:
    - permutations_results-50-500 (sizes 50-500)
    - permutation_results-550-1000 (sizes 550-1000)
    - permutation_results-FirstRun-50-to-1000 (all sizes 50-1000)
    
    Includes:
    - Complete coverage libraries (100%)
    - High partial coverage libraries (85-90%)
    Excludes medium coverage libraries (< 85%)
    """
    all_sizes = list(range(50, 1050, 50))
    library_sets = {}
    base_path = 'results'
    
    # Check all three folders and merge library sets
    for size in all_sizes:
        libraries_at_size = set()
        
        # Check 50-500 folder
        if size < 550:
            size_dir = f'{base_path}/permutations_results-50-500/size_{size}'
            if os.path.exists(size_dir):
                libraries_at_size.update(get_libraries_in_size(size_dir))
        
        # Check 550-1000 folder
        if size >= 550:
            size_dir = f'{base_path}/permutation_results-550-1000/size_{size}'
            if os.path.exists(size_dir):
                libraries_at_size.update(get_libraries_in_size(size_dir))
        
        # Check FirstRun folder (has all sizes)
        size_dir = f'{base_path}/permutation_results-FirstRun-50-to-1000/size_{size}'
        if os.path.exists(size_dir):
            libraries_at_size.update(get_libraries_in_size(size_dir))
        
        if libraries_at_size:
            library_sets[size] = libraries_at_size
    
    # Count coverage for each library
    library_coverage = defaultdict(set)
    
    for size, libraries in library_sets.items():
        for lib in libraries:
            library_coverage[lib].add(size)
    
    # Filter libraries by coverage threshold
    total_sizes = len(all_sizes)
    high_coverage_libraries = set()
    
    for library, sizes_present in library_coverage.items():
        coverage = len(sizes_present) / total_sizes
        if coverage >= min_coverage:
            high_coverage_libraries.add(library)
    
    return high_coverage_libraries, library_sets

def combine_library_data_for_size(size, library, output_dir):
    """Combine all permutation files for a specific library and size from all folders.
    
    Checks all three folders and merges data:
    - permutations_results-50-500 (sizes 50-500)
    - permutation_results-550-1000 (sizes 550-1000)
    - permutation_results-FirstRun-50-to-1000 (all sizes)
    
    Handles permutation numbering to ensure uniqueness across folders.
    """
    all_data = []
    permutation_count = 0
    max_perm_num = 0  # Track max permutation number to ensure uniqueness
    
    # List of folders to check (in order)
    folders_to_check = []
    
    # Check 50-500 folder
    if size < 550:
        size_dir = f'results/permutations_results-50-500/size_{size}'
        if os.path.exists(size_dir):
            folders_to_check.append(size_dir)
    
    # Check 550-1000 folder
    if size >= 550:
        size_dir = f'results/permutation_results-550-1000/size_{size}'
        if os.path.exists(size_dir):
            folders_to_check.append(size_dir)
    
    # Check FirstRun folder (has all sizes)
    size_dir = f'results/permutation_results-FirstRun-50-to-1000/size_{size}'
    if os.path.exists(size_dir):
        folders_to_check.append(size_dir)
    
    # Process each folder
    for size_dir in folders_to_check:
        # Get all permutation files
        permutation_files = sorted([f for f in os.listdir(size_dir) if f.endswith('.tsv')])
        
        for perm_file in permutation_files:
            file_path = os.path.join(size_dir, perm_file)
            try:
                df = pd.read_csv(file_path, sep='\t')
                # Filter for this library
                library_df = df[df['Library'] == library].copy()
                if not library_df.empty:
                    # Extract permutation number from filename
                    perm_num_str = perm_file.replace('permutation_', '').replace('.tsv', '')
                    perm_num = int(perm_num_str)
                    
                    # Ensure unique permutation numbers across folders
                    # Add offset based on folder to avoid collisions
                    if 'FirstRun' in size_dir:
                        perm_num_offset = 10000  # FirstRun gets 10000+
                    elif '550-1000' in size_dir:
                        perm_num_offset = 20000  # 550-1000 gets 20000+
                    else:
                        perm_num_offset = 0  # 50-500 gets 0-9999
                    
                    library_df['Permutation'] = perm_num + perm_num_offset
                    library_df['Source_Folder'] = os.path.basename(os.path.dirname(size_dir))
                    all_data.append(library_df)
                    permutation_count += 1
                    max_perm_num = max(max_perm_num, perm_num + perm_num_offset)
            except Exception as e:
                print(f"Warning: Error reading {file_path}: {e}")
                continue
    
    if all_data:
        # Combine all data
        combined_df = pd.concat(all_data, ignore_index=True)
        
        # Create output directory
        size_output_dir = os.path.join(output_dir, f'size_{size}')
        os.makedirs(size_output_dir, exist_ok=True)
        
        # Create safe filename from library name
        safe_library_name = library.replace('/', '_').replace('\\', '_').replace(':', '_')
        output_file = os.path.join(size_output_dir, f'{safe_library_name}.tsv')
        
        # Save combined file
        combined_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Combined {permutation_count} permutations for {library} at size {size}")
    
    return permutation_count

def main():
    """Main function to combine permutation data."""
    print("Finding libraries with high coverage (>= 85%)...")
    print("Including complete coverage (100%) and high partial coverage (85-90%) libraries.")
    print("Excluding medium coverage libraries (< 85%).")
    high_coverage_libraries, library_sets = find_libraries_with_high_coverage(min_coverage=0.85)
    
    print(f"\nFound {len(high_coverage_libraries)} libraries with >= 85% coverage:")
    for lib in sorted(high_coverage_libraries):
        sizes_count = len([s for s, libs in library_sets.items() if lib in libs])
        coverage_pct = (sizes_count / 20) * 100
        print(f"  - {lib} ({coverage_pct:.1f}% coverage, {sizes_count}/20 sizes)")
    
    # Create output directory
    output_base = 'results/combined_permutation_data'
    os.makedirs(output_base, exist_ok=True)
    
    # Process each size and library
    all_sizes = list(range(50, 1050, 50))
    summary_data = []
    
    print(f"\nCombining data for {len(high_coverage_libraries)} libraries across {len(all_sizes)} sizes...")
    
    for size in all_sizes:
        print(f"\nProcessing size {size}...")
        for library in sorted(high_coverage_libraries):
            # Only process if library has data at this size
            if size in library_sets and library in library_sets[size]:
                perm_count = combine_library_data_for_size(size, library, output_base)
                summary_data.append({
                    'Library': library,
                    'Gene_List_Size': size,
                    'Number_of_Permutations': perm_count
                })
            else:
                # Library doesn't have data at this size (missing due to no significant results)
                summary_data.append({
                    'Library': library,
                    'Gene_List_Size': size,
                    'Number_of_Permutations': 0
                })
    
    # Create summary table
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_base, 'permutation_summary_table.tsv')
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"\nSummary table saved to: {summary_file}")
    
    # Create pivot table for better readability
    pivot_df = summary_df.pivot_table(
        index='Library',
        columns='Gene_List_Size',
        values='Number_of_Permutations',
        fill_value=0
    )
    pivot_file = os.path.join(output_base, 'permutation_summary_pivot.tsv')
    pivot_df.to_csv(pivot_file, sep='\t')
    print(f"Pivot table saved to: {pivot_file}")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"\nTotal libraries with high coverage (>= 85%): {len(high_coverage_libraries)}")
    print(f"Total gene list sizes: {len(all_sizes)}")
    print(f"\nPermutations per library per size:")
    print(pivot_df.to_string())
    
    # Calculate coverage stats for metadata
    library_coverage_stats = {}
    for library in high_coverage_libraries:
        sizes_present = [s for s, libs in library_sets.items() if library in libs]
        library_coverage_stats[library] = {
            'sizes_present': sorted(sizes_present),
            'count': len(sizes_present),
            'coverage_percent': (len(sizes_present) / len(all_sizes)) * 100
        }
    
    # Save metadata
    metadata = {
        'libraries': sorted(list(high_coverage_libraries)),
        'library_coverage': library_coverage_stats,
        'gene_list_sizes': all_sizes,
        'total_libraries': len(high_coverage_libraries),
        'total_sizes': len(all_sizes),
        'min_coverage_threshold': 0.85,
        'output_directory': output_base
    }
    metadata_file = os.path.join(output_base, 'metadata.json')
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"\nMetadata saved to: {metadata_file}")

if __name__ == '__main__':
    main()

