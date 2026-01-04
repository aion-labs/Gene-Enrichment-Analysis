#!/bin/bash
# Script to rebuild parquet files for benchmarking with all 13 libraries

set -e

echo "=========================================="
echo "Rebuilding Parquet Files for Benchmarking"
echo "=========================================="
echo ""

# Step 1: Extract cluster statistics from combined permutation data
echo "Step 1: Extracting cluster statistics from combined permutation data..."
echo "This may take several minutes (processing ~40,000 permutations)..."
python extract_cluster_stats_from_combined_data.py \
    --combined-data-dir results/combined_permutation_data \
    --output results/permutation_cluster_statistics.tsv

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to extract cluster statistics"
    exit 1
fi

echo ""
echo "✓ Cluster statistics extracted successfully"
echo ""

# Step 2: Prepare parquet files
echo "Step 2: Preparing parquet files from cluster statistics..."
python prepare_parquet_cluster_stats.py \
    --input results/permutation_cluster_statistics.tsv \
    --output-dir results/permutation_cluster_statistics_parquet

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to prepare parquet files"
    exit 1
fi

echo ""
echo "✓ Parquet files prepared successfully"
echo ""
echo "=========================================="
echo "Summary:"
echo "=========================================="
echo "Cluster statistics: results/permutation_cluster_statistics.tsv"
echo "Parquet files: results/permutation_cluster_statistics_parquet/"
echo ""
echo "Libraries included (13 total):"
echo "  - C2: CGP: Chemical & genetic perturbations"
echo "  - C2: CP: BioCarta"
echo "  - C2: CP: Canonical pathways"
echo "  - C2: CP: KEGG MEDICUS"
echo "  - C2: CP: Pathway Interaction Database"
echo "  - C2: CP: WikiPathways"
echo "  - GO BP"
echo "  - GO CC"
echo "  - GO MF"
echo "  - H: Hallmark Gene Sets"
echo "  - KEGG"
echo "  - Protein Interaction"
echo "  - Reactome"
echo ""

