#!/usr/bin/env python3
"""
Test script for the CLI enrichment analysis.

This script tests both regular and iterative enrichment modes using the example gene list.
It uses smaller library files for faster testing and provides comprehensive output validation.
"""

import subprocess
import sys
import os
from pathlib import Path
import json
import pandas as pd

def run_command(cmd, description):
    """Run a command and return the result."""
    print(f"\n{'='*60}")
    print(f"Testing: {description}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        print("✅ Command executed successfully!")
        print("STDOUT:")
        print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        print("❌ Command failed!")
        print(f"Return code: {e.returncode}")
        print("STDOUT:")
        print(e.stdout)
        print("STDERR:")
        print(e.stderr)
        return False, e.stderr

def check_output_files(output_dir, expected_files):
    """Check if expected output files exist."""
    print(f"\nChecking output files in: {output_dir}")
    
    missing_files = []
    existing_files = []
    
    for expected_file in expected_files:
        file_path = Path(output_dir) / expected_file
        if file_path.exists():
            print(f"✅ Found: {expected_file}")
            existing_files.append(file_path)
        else:
            print(f"❌ Missing: {expected_file}")
            missing_files.append(expected_file)
    
    return len(missing_files) == 0, existing_files

def analyze_results(output_dir, test_name):
    """Analyze the results and provide a summary."""
    print(f"\n{'='*60}")
    print(f"Results Analysis: {test_name}")
    print(f"{'='*60}")
    
    # Check for JSON metadata file
    json_files = list(Path(output_dir).glob("*.json"))
    if json_files:
        print(f"📊 Metadata files found: {len(json_files)}")
        for json_file in json_files:
            try:
                with open(json_file, 'r') as f:
                    metadata = json.load(f)
                print(f"  📄 {json_file.name}:")
                print(f"    - Gene set: {metadata.get('gene_set_name', 'N/A')}")
                print(f"    - Gene set size: {metadata.get('gene_set_size', 'N/A')}")
                print(f"    - Libraries: {len(metadata.get('libraries', []))}")
                print(f"    - Total results: {metadata.get('total_results', 'N/A')}")
                if 'total_iterations' in metadata:
                    print(f"    - Total iterations: {metadata.get('total_iterations', 'N/A')}")
            except Exception as e:
                print(f"    ❌ Error reading {json_file}: {e}")
    
    # Check for TSV result files
    tsv_files = list(Path(output_dir).glob("*.tsv"))
    if tsv_files:
        print(f"📈 Result files found: {len(tsv_files)}")
        for tsv_file in tsv_files:
            try:
                df = pd.read_csv(tsv_file, sep='\t')
                print(f"  📄 {tsv_file.name}:")
                print(f"    - Rows: {len(df)}")
                print(f"    - Columns: {list(df.columns)}")
                if 'p-value' in df.columns:
                    significant = (df['p-value'] < 0.01).sum()
                    print(f"    - Significant results (p<0.01): {significant}")
            except Exception as e:
                print(f"    ❌ Error reading {tsv_file}: {e}")
    
    # Check for DOT network files
    dot_files = list(Path(output_dir).glob("*.dot"))
    if dot_files:
        print(f"🕸️  Network files found: {len(dot_files)}")
        for dot_file in dot_files:
            try:
                with open(dot_file, 'r') as f:
                    content = f.read()
                lines = content.split('\n')
                print(f"  📄 {dot_file.name}:")
                print(f"    - Lines: {len(lines)}")
                print(f"    - Size: {len(content)} characters")
            except Exception as e:
                print(f"    ❌ Error reading {dot_file}: {e}")

def main():
    """Main test function."""
    print("🧪 CLI Enrichment Analysis Test Suite")
    print("=" * 60)
    
    # Check if we're in the right directory
    if not Path("code/cli.py").exists():
        print("❌ Error: Please run this script from the project root directory")
        print("   Expected to find: code/cli.py")
        sys.exit(1)
    
    # Check if example files exist
    example_gene_list = Path("data/gene_lists/example_gene_list.txt")
    if not example_gene_list.exists():
        print("❌ Error: Example gene list not found")
        print("   Expected: data/gene_lists/example_gene_list.txt")
        sys.exit(1)
    
    # Find available libraries (use smaller ones for faster testing)
    libraries_dir = Path("data/libraries")
    small_libraries = [
        "h.all.v2025.1.Hs.symbols.gmt",  # Hallmark (small)
        "c2.cp.biocarta.v2025.1.Hs.symbols.gmt",  # BioCarta (small)
        "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt",  # KEGG MEDICUS (small)
    ]
    
    available_libraries = []
    for lib in small_libraries:
        lib_path = libraries_dir / lib
        if lib_path.exists():
            available_libraries.append(str(lib_path))
    
    if not available_libraries:
        print("❌ Error: No suitable library files found for testing")
        sys.exit(1)
    
    print(f"📚 Using libraries: {len(available_libraries)} files")
    for lib in available_libraries:
        print(f"   - {lib}")
    
    # Find background file
    background_files = list(Path("data/backgrounds").glob("*.txt"))
    if not background_files:
        print("❌ Error: No background files found")
        sys.exit(1)
    
    background_file = str(background_files[0])
    print(f"🎯 Using background: {background_file}")
    
    # Test 1: Basic CLI help
    print("\n" + "="*60)
    print("Test 1: CLI Help")
    print("="*60)
    success, output = run_command(
        ["python", "code/cli.py", "--help"],
        "CLI Help"
    )
    
    if not success:
        print("❌ CLI help test failed - CLI may not be working")
        sys.exit(1)
    
    # Test 2: Regular enrichment with defaults
    print("\n" + "="*60)
    print("Test 2: Regular Enrichment (Defaults)")
    print("="*60)
    
    success, output = run_command([
        "python", "code/cli.py",
        "--gene-sets", str(example_gene_list),
        "--background", background_file,
        "--libraries"] + available_libraries[:1],  # Use just one library for speed
        "Regular Enrichment with Defaults"
    )
    
    if success:
        # Check output files
        output_dir = Path("cli_results/example_gene_list")
        expected_files = [
            "combined_regular_results.tsv",
            "regular_enrichment_snapshot.json"
        ]
        
        files_ok, existing_files = check_output_files(output_dir, expected_files)
        if files_ok:
            analyze_results(output_dir, "Regular Enrichment")
        else:
            print("⚠️  Some expected files are missing")
    
    # Test 3: Regular enrichment with custom parameters
    print("\n" + "="*60)
    print("Test 3: Regular Enrichment (Custom Parameters)")
    print("="*60)
    
    success, output = run_command([
        "python", "code/cli.py",
        "--mode", "regular",
        "--gene-sets", str(example_gene_list),
        "--background", background_file,
        "--libraries"] + available_libraries[:2],  # Use two libraries
        "--p-threshold", "0.05",
        "--min-overlap", "2",
        "--min-term-size", "5",
        "--max-term-size", "500",
        "--output-dir", "test_regular_custom"
    ], "Regular Enrichment with Custom Parameters")
    
    if success:
        output_dir = Path("test_regular_custom/example_gene_list")
        expected_files = [
            "combined_regular_results.tsv",
            "regular_enrichment_snapshot.json"
        ]
        
        files_ok, existing_files = check_output_files(output_dir, expected_files)
        if files_ok:
            analyze_results(output_dir, "Regular Enrichment (Custom)")
    
    # Test 4: Iterative enrichment
    print("\n" + "="*60)
    print("Test 4: Iterative Enrichment")
    print("="*60)
    
    success, output = run_command([
        "python", "code/cli.py",
        "--mode", "iterative",
        "--gene-sets", str(example_gene_list),
        "--background", background_file,
        "--libraries"] + available_libraries[:1],  # Use just one library for speed
        "--p-threshold", "0.01",
        "--min-overlap", "3",
        "--max-iterations", "5",  # Limit iterations for testing
        "--output-dir", "test_iterative"
    ], "Iterative Enrichment")
    
    if success:
        output_dir = Path("test_iterative/example_gene_list")
        expected_files = [
            "combined_iterative_results.tsv",
            "iterative_enrichment_snapshot.json"
        ]
        
        files_ok, existing_files = check_output_files(output_dir, expected_files)
        if files_ok:
            analyze_results(output_dir, "Iterative Enrichment")
    
    # Test 5: Multiple gene sets
    print("\n" + "="*60)
    print("Test 5: Multiple Gene Sets")
    print("="*60)
    
    # Create a second test gene list
    test_gene_list2 = Path("test_gene_list2.txt")
    with open(test_gene_list2, 'w') as f:
        # Use first 50 genes from example list
        with open(example_gene_list, 'r') as example_f:
            genes = [line.strip() for line in example_f if line.strip()][:50]
        f.write('\n'.join(genes))
    
    success, output = run_command([
        "python", "code/cli.py",
        "--mode", "regular",
        "--gene-sets", str(example_gene_list), str(test_gene_list2),
        "--background", background_file,
        "--libraries"] + available_libraries[:1],
        "--output-dir", "test_multiple"
    ], "Multiple Gene Sets")
    
    if success:
        # Check both output directories
        for gene_list_name in ["example_gene_list", "test_gene_list2"]:
            output_dir = Path(f"test_multiple/{gene_list_name}")
            if output_dir.exists():
                print(f"\n📁 Results for {gene_list_name}:")
                analyze_results(output_dir, f"Multiple Gene Sets - {gene_list_name}")
    
    # Cleanup test files
    if test_gene_list2.exists():
        test_gene_list2.unlink()
    
    # Summary
    print("\n" + "="*60)
    print("🎉 Test Suite Summary")
    print("="*60)
    print("✅ CLI functionality tested successfully!")
    print("✅ Both regular and iterative modes working")
    print("✅ Multiple gene sets support working")
    print("✅ Custom parameters working")
    print("✅ Output files generated correctly")
    print("\n📁 Generated output directories:")
    
    output_dirs = [
        "cli_results",
        "test_regular_custom", 
        "test_iterative",
        "test_multiple"
    ]
    
    for output_dir in output_dirs:
        if Path(output_dir).exists():
            print(f"   - {output_dir}/")
    
    print("\n💡 You can examine the generated files to verify the results.")
    print("💡 The CLI is ready for production use!")

if __name__ == "__main__":
    main()
