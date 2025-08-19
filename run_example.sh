#!/bin/bash

# Example script for running the CLI enrichment analysis
# This script demonstrates various ways to use the CLI with the example gene list

set -e  # Exit on any error

echo "🧬 CLI Enrichment Analysis - Example Script"
echo "=========================================="

# Check if we're in the right directory
if [ ! -f "code/cli.py" ]; then
    echo "❌ Error: Please run this script from the project root directory"
    echo "   Expected to find: code/cli.py"
    exit 1
fi

# Check if virtual environment exists
if [ ! -d ".venv" ]; then
    echo "❌ Error: Virtual environment not found. Please run:"
    echo "   python -m venv .venv"
    echo "   source .venv/bin/activate"
    echo "   pip install -r requirements.txt"
    exit 1
fi

# Activate virtual environment
echo "🔧 Activating virtual environment..."
source .venv/bin/activate

# Check if example files exist
if [ ! -f "data/gene_lists/example_gene_list.txt" ]; then
    echo "❌ Error: Example gene list not found"
    echo "   Expected: data/gene_lists/example_gene_list.txt"
    exit 1
fi

echo "✅ Environment ready!"
echo ""

# Create output directory for examples
EXAMPLE_OUTPUT_DIR="example_results_$(date +%Y%m%d_%H%M%S)"
echo "📁 Output directory: $EXAMPLE_OUTPUT_DIR"
echo ""

# Example 1: Basic usage with defaults (all libraries)
echo "🔬 Example 1: Basic Analysis (All Libraries)"
echo "--------------------------------------------"
echo "Running: python code/cli.py --gene-sets data/gene_lists/example_gene_list.txt"
echo ""

python code/cli.py \
    --gene-sets data/gene_lists/example_gene_list.txt \
    --output-dir "${EXAMPLE_OUTPUT_DIR}/example1_basic"

echo ""
echo "✅ Example 1 completed!"
echo ""

# Example 2: Regular enrichment with specific libraries
echo "🔬 Example 2: Regular Enrichment (Specific Libraries)"
echo "----------------------------------------------------"
echo "Running: python code/cli.py --mode regular --gene-sets data/gene_lists/example_gene_list.txt --libraries h.all.v2025.1.Hs.symbols.gmt c5.go.bp.v2025.1.Hs.symbols.gmt"
echo ""

python code/cli.py \
    --mode regular \
    --gene-sets data/gene_lists/example_gene_list.txt \
    --libraries data/libraries/h.all.v2025.1.Hs.symbols.gmt data/libraries/c5.go.bp.v2025.1.Hs.symbols.gmt \
    --p-threshold 0.01 \
    --min-overlap 3 \
    --output-dir "${EXAMPLE_OUTPUT_DIR}/example2_regular"

echo ""
echo "✅ Example 2 completed!"
echo ""

# Example 3: Iterative enrichment with limited iterations
echo "🔬 Example 3: Iterative Enrichment (Limited Iterations)"
echo "------------------------------------------------------"
echo "Running: python code/cli.py --mode iterative --gene-sets data/gene_lists/example_gene_list.txt --max-iterations 5"
echo ""

python code/cli.py \
    --mode iterative \
    --gene-sets data/gene_lists/example_gene_list.txt \
    --libraries data/libraries/h.all.v2025.1.Hs.symbols.gmt \
    --p-threshold 0.01 \
    --min-overlap 3 \
    --max-iterations 5 \
    --output-dir "${EXAMPLE_OUTPUT_DIR}/example3_iterative"

echo ""
echo "✅ Example 3 completed!"
echo ""

# Example 4: Custom parameters with different p-value method
echo "🔬 Example 4: Custom Parameters (Hypergeometric Test)"
echo "----------------------------------------------------"
echo "Running: python code/cli.py --mode regular --gene-sets data/gene_lists/example_gene_list.txt --method 'Hypergeometric Test' --p-threshold 0.05"
echo ""

python code/cli.py \
    --mode regular \
    --gene-sets data/gene_lists/example_gene_list.txt \
    --libraries data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt \
    --method "Hypergeometric Test" \
    --p-threshold 0.05 \
    --min-overlap 2 \
    --min-term-size 10 \
    --max-term-size 500 \
    --output-dir "${EXAMPLE_OUTPUT_DIR}/example4_custom"

echo ""
echo "✅ Example 4 completed!"
echo ""

# Example 5: Multiple gene sets
echo "🔬 Example 5: Multiple Gene Sets"
echo "--------------------------------"
echo "Creating a second test gene list and running analysis on both..."
echo ""

# Create a second test gene list (first 50 genes from example)
head -50 data/gene_lists/example_gene_list.txt > temp_gene_list2.txt

python code/cli.py \
    --mode regular \
    --gene-sets data/gene_lists/example_gene_list.txt temp_gene_list2.txt \
    --libraries data/libraries/h.all.v2025.1.Hs.symbols.gmt \
    --p-threshold 0.01 \
    --output-dir "${EXAMPLE_OUTPUT_DIR}/example5_multiple"

# Clean up temporary file
rm -f temp_gene_list2.txt

echo ""
echo "✅ Example 5 completed!"
echo ""

# Summary
echo "🎉 All Examples Completed Successfully!"
echo "======================================="
echo ""
echo "📁 Generated output directories:"
echo "   - ${EXAMPLE_OUTPUT_DIR}/example1_basic/     (Basic analysis with all libraries)"
echo "   - ${EXAMPLE_OUTPUT_DIR}/example2_regular/   (Regular mode with specific libraries)"
echo "   - ${EXAMPLE_OUTPUT_DIR}/example3_iterative/ (Iterative mode with limited iterations)"
echo "   - ${EXAMPLE_OUTPUT_DIR}/example4_custom/    (Custom parameters and method)"
echo "   - ${EXAMPLE_OUTPUT_DIR}/example5_multiple/  (Multiple gene sets)"
echo ""
echo "📊 Each directory contains:"
echo "   - Individual library result files (*_regular_results.tsv or *_iterative_results.tsv)"
echo "   - Combined results file (combined_*_results.tsv)"
echo "   - Metadata file (*_enrichment_snapshot.json)"
echo "   - Network files (*_network.dot) for iterative mode"
echo ""
echo "🔍 To examine results:"
echo "   ls -la ${EXAMPLE_OUTPUT_DIR}/"
echo "   head -20 ${EXAMPLE_OUTPUT_DIR}/example1_basic/example_gene_list/combined_regular_results.tsv"
echo ""
echo "💡 CLI Usage Tips:"
echo "   - Use --help to see all options: python code/cli.py --help"
echo "   - Omit --libraries to use all available libraries"
echo "   - Use --output-dir to specify custom output location"
echo "   - Adjust --p-threshold and --min-overlap for different sensitivity"
echo ""
echo "✨ The CLI is ready for your own gene lists!"
