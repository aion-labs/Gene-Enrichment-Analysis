#!/usr/bin/env python3
"""
Analyze which libraries have permutation data at which gene list sizes.
"""

import pandas as pd
import os
from collections import defaultdict
from pathlib import Path

def get_libraries_in_size(size_dir):
    """Get unique libraries from a size directory by checking all permutation files."""
    libraries = set()
    
    if not os.path.exists(size_dir):
        return libraries
    
    # Get all permutation files
    perm_files = [f for f in os.listdir(size_dir) if f.endswith('.tsv')]
    
    if not perm_files:
        return libraries
    
    # Check first few files to get all libraries (they should all have the same libraries)
    # But to be safe, check multiple files
    files_to_check = perm_files[:min(10, len(perm_files))]
    
    for perm_file in files_to_check:
        file_path = os.path.join(size_dir, perm_file)
        try:
            df = pd.read_csv(file_path, sep='\t')
            libraries.update(df['Library'].unique())
        except Exception as e:
            print(f"Warning: Error reading {file_path}: {e}")
            continue
    
    return libraries

def main():
    """Analyze library coverage across all sizes."""
    all_sizes = list(range(50, 1050, 50))
    
    # Dictionary: library -> set of sizes where it exists
    library_coverage = defaultdict(set)
    
    # Dictionary: size -> set of libraries
    size_libraries = {}
    
    print("Scanning permutation folders...")
    print("="*80)
    
    # Check 50-500 folder
    for size in range(50, 550, 50):
        size_dir = f'results/permutations_results-50-500/size_{size}'
        libraries = get_libraries_in_size(size_dir)
        size_libraries[size] = libraries
        for lib in libraries:
            library_coverage[lib].add(size)
        print(f"Size {size:4d}: {len(libraries):2d} libraries")
    
    # Check 550-1000 folder
    for size in range(550, 1050, 50):
        size_dir = f'results/permutation_results-550-1000/size_{size}'
        libraries = get_libraries_in_size(size_dir)
        size_libraries[size] = libraries
        for lib in libraries:
            library_coverage[lib].add(size)
        print(f"Size {size:4d}: {len(libraries):2d} libraries")
    
    print("="*80)
    print(f"\nTotal unique libraries found: {len(library_coverage)}")
    print("\n" + "="*80)
    print("LIBRARY COVERAGE REPORT")
    print("="*80)
    
    # Sort libraries by name
    sorted_libraries = sorted(library_coverage.keys())
    
    # Create detailed report
    report_lines = []
    report_lines.append(f"{'Library':<60} {'Sizes Present':<20} {'Count':<10} {'Coverage'}")
    report_lines.append("-" * 120)
    
    for library in sorted_libraries:
        sizes = sorted(library_coverage[library])
        size_count = len(sizes)
        size_range = f"{min(sizes)}-{max(sizes)}" if sizes else "None"
        
        # Calculate coverage percentage (out of 20 possible sizes)
        coverage_pct = (size_count / 20) * 100
        
        # Create size list string (abbreviated if too long)
        if size_count <= 10:
            sizes_str = ",".join(map(str, sizes))
        else:
            sizes_str = f"{sizes[0]}-{sizes[-1]} ({size_count} sizes)"
        
        report_lines.append(f"{library:<60} {sizes_str:<20} {size_count:<10} {coverage_pct:.1f}%")
    
    print("\n".join(report_lines))
    
    # Create a matrix showing which libraries have data at which sizes
    print("\n" + "="*80)
    print("COVERAGE MATRIX (X = data exists)")
    print("="*80)
    
    # Create matrix header
    matrix_lines = []
    header = f"{'Library':<60} " + " ".join([f"{s:>4}" for s in all_sizes])
    matrix_lines.append(header)
    matrix_lines.append("-" * len(header))
    
    for library in sorted_libraries:
        row = f"{library:<60} "
        for size in all_sizes:
            if size in library_coverage[library]:
                row += "   X"
            else:
                row += "   ."
        matrix_lines.append(row)
    
    print("\n".join(matrix_lines))
    
    # Save detailed report to file
    output_file = 'results/library_coverage_report.txt'
    os.makedirs('results', exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("LIBRARY COVERAGE ANALYSIS\n")
        f.write("="*80 + "\n\n")
        f.write("This report shows which libraries have permutation data at which gene list sizes.\n\n")
        f.write("\n".join(report_lines))
        f.write("\n\n")
        f.write("="*80 + "\n")
        f.write("COVERAGE MATRIX (X = data exists, . = missing)\n")
        f.write("="*80 + "\n\n")
        f.write("\n".join(matrix_lines))
        f.write("\n\n")
        f.write("="*80 + "\n")
        f.write("DETAILED SIZE INFORMATION\n")
        f.write("="*80 + "\n\n")
        
        for library in sorted_libraries:
            sizes = sorted(library_coverage[library])
            f.write(f"{library}:\n")
            f.write(f"  Sizes: {', '.join(map(str, sizes))}\n")
            f.write(f"  Count: {len(sizes)}/20 sizes\n")
            f.write(f"  Coverage: {(len(sizes)/20)*100:.1f}%\n")
            f.write(f"  Range: {min(sizes)} to {max(sizes)}\n")
            
            # Identify gaps
            all_sizes_set = set(all_sizes)
            missing = sorted(all_sizes_set - library_coverage[library])
            if missing:
                f.write(f"  Missing sizes: {', '.join(map(str, missing))}\n")
            else:
                f.write(f"  Complete coverage (all sizes present)\n")
            f.write("\n")
    
    print(f"\nDetailed report saved to: {output_file}")
    
    # Create CSV summary
    csv_data = []
    for library in sorted_libraries:
        sizes = sorted(library_coverage[library])
        csv_data.append({
            'Library': library,
            'Sizes_Present': ','.join(map(str, sizes)),
            'Size_Count': len(sizes),
            'Coverage_Percent': (len(sizes)/20)*100,
            'Min_Size': min(sizes) if sizes else None,
            'Max_Size': max(sizes) if sizes else None,
            'Complete': len(sizes) == 20
        })
    
    df = pd.DataFrame(csv_data)
    csv_file = 'results/library_coverage_summary.csv'
    df.to_csv(csv_file, index=False)
    print(f"CSV summary saved to: {csv_file}")

if __name__ == '__main__':
    main()

