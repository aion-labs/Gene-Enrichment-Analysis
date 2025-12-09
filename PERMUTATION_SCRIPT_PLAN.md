# Permutation Testing Script Plan

## Overview
Generate a distribution of p-values from random gene lists using iGEA (iterative Gene Enrichment Analysis) for statistical validation in a research paper. This script will run offline and generate 20,000 TSV files (20 gene list sizes × 1000 permutations each).

## Libraries
The following 5 libraries will be used:
1. **Reactome**: `c2.cp.reactome.v2025.1.Hs.symbols.gmt`
2. **KEGG**: `c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt`
3. **GO Biological Process**: `c5.go.bp.v2025.1.Hs.symbols.gmt`
4. **GO Molecular Function**: `c5.go.mf.v2025.1.Hs.symbols.gmt`
5. **GO Cellular Component**: `c5.go.cc.v2025.1.Hs.symbols.gmt`

## Parameters
- **Background**: `all_genes.txt` from `data/backgrounds/`
- **P-value threshold**: 0.05 (raw p-value, as specified)
- **Min overlap**: 3 (default from Streamlit)
- **Min term size**: 10 (default from Streamlit)
- **Max term size**: 600 (default from Streamlit)
- **Max iterations**: 10 (default from Streamlit)
- **P-value method**: "Fisher's Exact Test" (default)
- **Gene format**: Gene Symbols

## Gene List Sizes
- Range: 50 to 1000 genes
- Increment: 50 genes
- Total sizes: 20 (50, 100, 150, ..., 1000)
- Permutations per size: 1000
- **Total runs**: 20,000 iGEA analyses

## Output Structure
```
permutation_results/
├── size_50/
│   ├── permutation_0001.tsv
│   ├── permutation_0002.tsv
│   ├── ...
│   └── permutation_1000.tsv
├── size_100/
│   ├── permutation_0001.tsv
│   ├── ...
│   └── permutation_1000.tsv
...
└── size_1000/
    ├── permutation_0001.tsv
    ├── ...
    └── permutation_1000.tsv
```

Each TSV file will contain:
- Combined results from all 5 libraries
- Same format as `to_dataframe()` output
- Columns: Library, Iteration, Term, Description, iteration overlapping genes, iteration p-value, iteration -log(p-value), Genes removed for next iteration, etc.

## Optimization Strategies

1. **Parallel Processing**:
   - Use `multiprocessing` to run multiple permutations simultaneously
   - Process multiple gene list sizes in parallel
   - Recommended: Use number of CPU cores available

2. **Batch Processing**:
   - Process libraries sequentially within each permutation (as they share the same gene set)
   - Group permutations into batches for better resource management

3. **Memory Management**:
   - Load libraries once and reuse across all permutations
   - Load background gene set once
   - Clear intermediate results after saving

4. **Progress Tracking**:
   - Log progress to file for resumability
   - Display progress bar/percentage
   - Save checkpoint after each completed permutation

5. **Error Handling**:
   - Continue on individual permutation failures
   - Log errors to separate file
   - Track success/failure counts

6. **Output Optimization**:
   - Write TSV files directly (no intermediate storage)
   - Use efficient TSV writing (pandas `to_csv` with tab separator)

## Implementation Steps

1. **Setup**:
   - Load background gene set (`all_genes.txt`)
   - Load all 5 libraries
   - Create output directory structure
   - Set up logging

2. **For each gene list size (50, 100, ..., 1000)**:
   - Create size-specific output directory
   - For each permutation (1 to 1000):
     - Randomly sample genes from background
     - Create GeneSet object
     - Run iGEA for each library sequentially
     - Combine results from all libraries
     - Save combined CSV file
     - Log progress

3. **Post-processing**:
   - Generate summary statistics
   - Create summary report of completed runs

## Expected Runtime
- Assuming ~10-30 seconds per permutation (5 libraries × 2-6 seconds each)
- 1000 permutations × 20 seconds average = ~20,000 seconds = ~5.5 hours per size
- 20 sizes × 5.5 hours = ~110 hours total (single-threaded)
- With 8-core parallelization: ~14 hours total
- With 16-core parallelization: ~7 hours total

## Resumability
- Check for existing TSV files before running permutation
- Skip already completed permutations
- Allow script to resume from last completed permutation

## Output TSV Format
Each TSV file will contain combined results from all libraries, with columns:
- Library
- Iteration
- Term
- Description
- iteration overlapping genes
- iteration p-value
- iteration -log(p-value)
- Genes removed for next iteration
- Full list overlapping genes (if available)
- Full list p-value (if available)
- Regular FDR (if available)
- Full list overlapping genes (gene names) (if available)

