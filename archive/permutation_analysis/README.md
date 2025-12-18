# Permutation Analysis Archive

This folder contains all files related to the permutation-based p-value correction analysis for iGEA.

## Structure

```
archive/permutation_analysis/
├── README.md (this file)
├── scripts/
│   ├── generate_permutation_distribution.py - Main permutation generation script
│   ├── code/generate_permutation_distribution.py - Code version of permutation script
│   ├── check_permutation_progress.py - Script to check permutation progress
│   ├── merge_permutation_results.py - Script to merge permutation results
│   └── pvalue_corrector.py - P-value correction implementation
├── results/
│   ├── permutation_results/ - Permutation data (1.25M rows)
│   │   ├── merged_permutation_results.tsv - Merged permutation results
│   │   └── Analyze Permutation Results.twbx - Tableau workbook
│   └── permutation_distribution.log - Log file from permutation generation
├── documentation/
│   ├── PERMUTATION_SCRIPT_PLAN.md - Initial plan for permutation script
│   ├── PERMUTATION_SCRIPT_README.md - README for permutation script
│   ├── P_VALUE_CORRECTION_RECOMMENDATIONS.md - Detailed recommendations for p-value correction
│   ├── P_VALUE_CORRECTION_SUMMARY.md - Summary of p-value correction approach
│   ├── PVALUE_CORRECTOR_INTEGRATION.md - Integration guide for p-value corrector
│   ├── GENE_LIST_SIZE_HANDLING.md - How gene list sizes are handled
│   └── STATISTICAL_CONCEPTS_EXPLAINED.md - Explanation of stratum and CDF concepts
└── simulation/
    ├── simulate_hiv_pvalue_correction.py - HIV gene list simulation script
    ├── example_pvalue_correction.py - Example usage script
    └── hiv_simulation/ - Results from HIV simulation
        ├── hiv_igea_with_correction.tsv - Results with corrected p-values
        └── hiv_igea_raw.tsv - Raw results for comparison
```

## Overview

This archive contains the complete permutation-based p-value correction system for iGEA:

1. **Permutation Generation**: Scripts to generate null distributions from random gene lists
2. **P-Value Correction**: Implementation of stratified empirical null distribution correction
3. **Documentation**: Comprehensive guides and explanations
4. **Simulation**: Example application on HIV gene list

## Key Concepts

- **Stratum**: A subgroup of data with the same characteristics (library, iteration, term size, overlap, gene list size)
- **CDF (Cumulative Distribution Function)**: Function that tells you what percentage of values are ≤ a given value
- **Stratified Correction**: Creating separate null distributions for each combination of factors

## Permutation Data

- **Size**: 1.25 million permutation results
- **Libraries**: GO BP, GO CC, GO MF, KEGG, Reactome
- **Gene List Sizes**: 50, 100, 150, ..., 1000 (increments of 50)
- **Iterations**: Up to 30 iterations per permutation
- **Strata**: 2,753 stratum CDFs built from permutation data

## Status

⚠️ **Note**: The p-value correction results showed unexpected behavior and require further investigation. This archive preserves all work for future reference and debugging.

## Date Archived

December 13, 2025
