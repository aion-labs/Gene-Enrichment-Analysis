#!/usr/bin/env python3
"""
Test script for running HIV dependency factors gene list through CLI.

This script tests the CLI with:
- HIV_dependency_factors.symbols.txt gene list
- All 11 libraries that have permutation data available
- P-value threshold: 0.01
- Default parameters for everything else
- Iterative mode with benchmarking enabled
"""

import subprocess
import sys
from pathlib import Path
from datetime import datetime

# Get project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# Test configuration
HIV_GENE_LIST = PROJECT_ROOT / "data" / "gene_lists" / "HIV_dependency_factors.symbols.txt"
BACKGROUND = PROJECT_ROOT / "data" / "backgrounds" / "all_genes.txt"
OUTPUT_DIR = Path(__file__).resolve().parent / "results"

# All 11 libraries with permutation data
LIBRARIES = [
    "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt",  # Reactome
    "data/libraries/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",  # KEGG
    "data/libraries/c5.go.bp.v2025.1.Hs.symbols.gmt",  # GO BP
    "data/libraries/c5.go.mf.v2025.1.Hs.symbols.gmt",  # GO MF
    "data/libraries/c5.go.cc.v2025.1.Hs.symbols.gmt",  # GO CC
    "data/libraries/c2.cp.biocarta.v2025.1.Hs.symbols.gmt",  # BioCarta
    "data/libraries/c2.cp.v2025.1.Hs.symbols.gmt",  # Canonical pathways
    "data/libraries/c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt",  # KEGG Medicus
    "data/libraries/c2.cp.pid.v2025.1.Hs.symbols.gmt",  # Pathway Interaction Database
    "data/libraries/c2.cp.wikipathways.v2025.1.Hs.symbols.gmt",  # WikiPathways
    "data/libraries/h.all.v2025.1.Hs.symbols.gmt",  # Hallmark
]

def check_prerequisites():
    """Check that all required files exist."""
    print("=" * 80)
    print("Checking Prerequisites")
    print("=" * 80)
    
    issues = []
    
    # Check gene list
    if not HIV_GENE_LIST.exists():
        issues.append(f"❌ HIV gene list not found: {HIV_GENE_LIST}")
    else:
        print(f"✅ HIV gene list found: {HIV_GENE_LIST}")
        with open(HIV_GENE_LIST, 'r') as f:
            num_genes = len([line for line in f if line.strip()])
        print(f"   Contains {num_genes} genes")
    
    # Check background
    if not BACKGROUND.exists():
        issues.append(f"❌ Background file not found: {BACKGROUND}")
    else:
        print(f"✅ Background file found: {BACKGROUND}")
    
    # Check libraries
    missing_libs = []
    for lib in LIBRARIES:
        lib_path = PROJECT_ROOT / lib
        if not lib_path.exists():
            missing_libs.append(lib)
    
    if missing_libs:
        issues.append(f"❌ Missing libraries: {missing_libs}")
    else:
        print(f"✅ All {len(LIBRARIES)} libraries found")
    
    # Check CLI
    cli_path = PROJECT_ROOT / "code" / "cli.py"
    if not cli_path.exists():
        issues.append(f"❌ CLI not found: {cli_path}")
    else:
        print(f"✅ CLI found: {cli_path}")
    
    if issues:
        print("\n❌ Prerequisites check failed:")
        for issue in issues:
            print(f"   {issue}")
        return False
    
    print("\n✅ All prerequisites met!")
    return True

def run_hiv_test():
    """Run the HIV gene list test."""
    print("\n" + "=" * 80)
    print("Running HIV Dependency Factors Test")
    print("=" * 80)
    print(f"\nConfiguration:")
    print(f"  Gene List: {HIV_GENE_LIST.name}")
    print(f"  Background: {BACKGROUND.name}")
    print(f"  Libraries: {len(LIBRARIES)} libraries with permutation data")
    print(f"  P-value threshold: 0.01")
    print(f"  Max iterations: 30")
    print(f"  Mode: iterative (with benchmarking)")
    print(f"  Output: {OUTPUT_DIR}")
    print()
    
    # Build command
    cmd = [
        sys.executable,
        str(PROJECT_ROOT / "code" / "cli.py"),
        "--genelist", str(HIV_GENE_LIST),
        "--background", str(BACKGROUND),
        "--mode", "iterative",
        "--p-threshold", "0.01",
        "--max-iterations", "30",  # Set max iterations to 30
        "--benchmark",  # Enable statistical benchmarking
        "--output-dir", str(OUTPUT_DIR),
    ]
    
    # Add all libraries (typer List[Path] requires --libraries before each library)
    for lib in LIBRARIES:
        cmd.extend(["--libraries", str(PROJECT_ROOT / lib)])
    
    print("Command:")
    print(" ".join(cmd))
    print("\n" + "-" * 80)
    print("Running analysis (this may take several minutes)...")
    print("-" * 80 + "\n")
    
    start_time = datetime.now()
    
    try:
        result = subprocess.run(
            cmd,
            cwd=str(PROJECT_ROOT),
            capture_output=False,  # Show output in real-time
            text=True,
            check=True
        )
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print("\n" + "=" * 80)
        print("✅ Test Completed Successfully!")
        print("=" * 80)
        print(f"Duration: {duration:.1f} seconds ({duration/60:.1f} minutes)")
        print(f"Results saved to: {OUTPUT_DIR}")
        
        return True
        
    except subprocess.CalledProcessError as e:
        print("\n" + "=" * 80)
        print("❌ Test Failed!")
        print("=" * 80)
        print(f"Return code: {e.returncode}")
        return False
    except KeyboardInterrupt:
        print("\n\n⚠️  Test interrupted by user")
        return False

def check_results():
    """Check that expected output files were generated."""
    print("\n" + "=" * 80)
    print("Checking Results")
    print("=" * 80)
    
    if not OUTPUT_DIR.exists():
        print(f"❌ Output directory not found: {OUTPUT_DIR}")
        return False
    
    # Find the gene set output directory (it will have a timestamp)
    gene_set_dirs = list(OUTPUT_DIR.glob("HIV_dependency_factors*"))
    if not gene_set_dirs:
        print(f"⚠️  No gene set output directory found in {OUTPUT_DIR}")
        return False
    
    results_dir = gene_set_dirs[0]
    print(f"✅ Results directory: {results_dir}")
    
    # Check for expected files
    expected_files = [
        "combined_iterative_results.tsv",
        "combined_network.dot",
    ]
    
    # Check for library-specific results
    library_results = list(results_dir.glob("*_iterative_results.tsv"))
    print(f"\n📊 Found {len(library_results)} library result files")
    
    # Check for benchmark files (if benchmarking was enabled)
    benchmark_files = [
        "statistical_benchmarks_table.tsv",
    ]
    # Also check for the statistical report (name depends on gene set name)
    report_files = list(results_dir.glob("*_statistical_report.txt"))
    
    found_files = []
    missing_files = []
    
    for expected in expected_files:
        file_path = results_dir / expected
        if file_path.exists():
            found_files.append(expected)
            size = file_path.stat().st_size
            print(f"  ✅ {expected} ({size:,} bytes)")
        else:
            missing_files.append(expected)
            print(f"  ⚠️  {expected} (not found)")
    
    # Check benchmark files
    for expected in benchmark_files:
        file_path = results_dir / expected
        if file_path.exists():
            found_files.append(expected)
            size = file_path.stat().st_size
            print(f"  ✅ {expected} ({size:,} bytes)")
        else:
            print(f"  ⚠️  {expected} (not found - benchmarking may have been skipped)")
    
    # Check for statistical report
    if report_files:
        for report_file in report_files:
            found_files.append(report_file.name)
            size = report_file.stat().st_size
            print(f"  ✅ {report_file.name} ({size:,} bytes)")
    else:
        print(f"  ⚠️  *_statistical_report.txt (not found - benchmarking may have been skipped)")
    
    if missing_files:
        print(f"\n⚠️  Some expected files are missing: {missing_files}")
        return False
    
    print(f"\n✅ All expected files found!")
    return True

def main():
    """Main test function."""
    print("\n" + "🧪 " * 30)
    print("HIV Dependency Factors - CLI Test")
    print("🧪 " * 30 + "\n")
    
    # Check prerequisites
    if not check_prerequisites():
        print("\n❌ Prerequisites check failed. Please fix issues above.")
        return 1
    
    # Run test
    success = run_hiv_test()
    
    if not success:
        return 1
    
    # Check results
    if not check_results():
        print("\n⚠️  Test completed but some expected files are missing.")
        return 1
    
    print("\n" + "=" * 80)
    print("🎉 ALL TESTS PASSED!")
    print("=" * 80)
    print(f"\nResults are available in: {OUTPUT_DIR}")
    print("\nYou can examine the results files to verify the analysis worked correctly.")
    print()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

