# Permutation Statistics Report Generator

## Overview

This script generates comprehensive statistical reports for permutation data across different gene list sizes and p-value thresholds.

## What It Does

1. **Loads merged permutation results** from the archive
2. **Processes each combination** of:
   - Gene list sizes: 50, 100, 150, 200, 250, 300, 400, 500, 750, 1000
   - P-value thresholds: 0.05, 0.01, 0.001
3. **For each permutation**:
   - Filters data by gene list size and p-value threshold
   - Builds network from enrichment results
   - Identifies clusters (connected components)
   - Extracts metrics for all clusters
4. **Computes statistics** (mean, std, min, max, median) for largest cluster metrics:
   - `largest_cluster_genes`: Number of genes in largest cluster
   - `largest_cluster_terms`: Number of terms in largest cluster
   - `largest_cluster_edges`: Number of edges in largest cluster
   - `largest_cluster_density`: Density of largest cluster
   - `largest_component_size`: Total size (genes + terms) of largest cluster
   - `fraction_in_largest_cluster`: Fraction of network in largest cluster
   - `fraction_edges_in_largest_cluster`: Fraction of edges in largest cluster
   - `largest_cluster_libraries`: Number of libraries in largest cluster
5. **Generates three output files**:
   - Summary report (text file)
   - Summary CSV (statistics by size/p-value)
   - Detailed CSV (all clusters from all permutations)

## Usage

```bash
python code/generate_permutation_statistics_report.py
```

## Output Files

All files are saved to `results/permutation_statistics/`:

1. **`permutation_statistics_report.txt`**
   - Human-readable summary report
   - Organized by gene list size and p-value threshold
   - Shows mean, std, min, max, median, and N for each metric

2. **`permutation_statistics_summary.csv`**
   - CSV format with one row per size/p-value combination
   - Columns: `gene_list_size`, `p_value_threshold`, `n_permutations`, and statistics columns
   - Statistics columns: `{metric}_mean`, `{metric}_std`, `{metric}_min`, `{metric}_max`, `{metric}_median`

3. **`permutation_clusters_detailed.csv`**
   - One row per cluster per permutation
   - Columns include:
     - `gene_list_size`: Gene list size for this permutation
     - `p_value_threshold`: P-value threshold used
     - `permutation_id`: Unique identifier for the permutation
     - `filename`: Original filename
     - `cluster_number`: Cluster number (1 = largest, 2 = second largest, etc.)
     - `n_genes`, `n_terms`, `n_edges`: Cluster-specific metrics
     - `cluster_size`: Total size (genes + terms)
     - `density`: Cluster density
     - `n_libraries`: Number of libraries in cluster
     - `largest_cluster_*`: Metrics for largest cluster (only populated for cluster_number == 1)
     - `n_connected_components`, `n_genes_total`, `n_terms_total`, `n_edges_total`: Network-wide metrics (only for largest cluster)

## Notes

- The script processes all permutations for each size/p-value combination
- Statistics are computed only for the **largest cluster** from each permutation
- The detailed CSV includes **all clusters** from all permutations
- Processing may take several minutes depending on the number of permutations
- Progress is logged to the console

## Example Output

### Summary Report Format
```
================================================================================
Gene List Size: 200
================================================================================

P-value threshold: 0.05
  Number of permutations: 1000

  largest_cluster_genes:
    Mean:   64.7234
    Std:    20.8123
    Min:    10.0000
    Max:    114.0000
    Median: 67.0000
    N:      1000
```

### Summary CSV Format
```csv
gene_list_size,p_value_threshold,n_permutations,largest_cluster_genes_mean,largest_cluster_genes_std,...
200,0.05,1000,64.7234,20.8123,...
200,0.01,850,58.2341,18.9234,...
```

### Detailed CSV Format
```csv
gene_list_size,p_value_threshold,permutation_id,filename,cluster_number,n_genes,n_terms,...
200,0.05,permutation_001,permutation_001.tsv,1,45,23,...
200,0.05,permutation_001,permutation_001.tsv,2,12,5,...
```

