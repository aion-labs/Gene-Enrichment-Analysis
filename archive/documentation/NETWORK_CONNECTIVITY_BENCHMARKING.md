# Network Connectivity Benchmarking for iGEA

## Overview

This document describes the network connectivity benchmarking approach for iGEA. The central hypothesis is that **real gene lists will have better network connectivity than random gene lists**, as measured by various graph metrics.

## Approach

### 1. Network Structure

iGEA results form a **bipartite graph**:
- **Gene nodes**: Individual genes from the input gene list
- **Term nodes**: Enriched biological terms/pathways
- **Edges**: Gene-term connections (gene is a member of that term)

### 2. Connectivity Metrics

The following metrics are computed to measure network connectivity:

#### Basic Metrics
- **n_genes**: Number of unique genes in the network
- **n_terms**: Number of enriched terms
- **n_edges**: Total number of gene-term connections

#### Connectivity Metrics
- **avg_connections_per_gene**: Average number of terms each gene is connected to
- **avg_connections_per_term**: Average number of genes in each term
- **network_density**: Ratio of actual edges to possible edges (n_edges / (n_genes × n_terms))

#### Graph Structure Metrics
- **n_connected_components**: Number of disconnected subgraphs
- **largest_component_size**: Size of the largest connected component
- **clustering_coefficient**: Measure of how interconnected the network is

#### Centrality Metrics
- **gene_centrality_mean/std/max**: Distribution of gene degrees (how many terms each gene connects to)
- **term_centrality_mean/std/max**: Distribution of term degrees (how many genes in each term)

### 3. Null Distribution from Permutations

Using the permutation results (random gene lists), we:

1. **Compute connectivity metrics** for each permutation
2. **Stratify by gene list size** (50, 100, 150, ..., 1000)
3. **Build null distributions** with statistics (mean, std, percentiles) for each metric

### 4. Benchmarking Real Results

For a real gene list:
1. Run iGEA to get results
2. Compute connectivity metrics
3. Compare against null distribution for the same gene list size
4. Calculate:
   - **Z-scores**: How many standard deviations above/below the null mean
   - **Percentiles**: What percentile of the null distribution the real value falls in
   - **Significance**: Whether connectivity is significantly better than random

## Usage

### Step 1: Compute Connectivity Metrics from Permutations

```bash
python code/network_connectivity_benchmark.py \
    --permutation-file archive/permutation_analysis/results/permutation_results/merged_permutation_results.tsv \
    --output-metrics results/connectivity_metrics.tsv \
    --output-null results/connectivity_null_distribution.json
```

This will:
- Process all permutation results
- Compute connectivity metrics for each permutation
- Build null distributions stratified by gene list size
- Save results to TSV and JSON files

### Step 2: Benchmark Real Results

```python
from code.network_connectivity_benchmark import (
    NetworkConnectivityAnalyzer,
    benchmark_real_results
)
import json

# Load null distribution
with open('results/connectivity_null_distribution.json', 'r') as f:
    null_distribution = json.load(f)

# Get iGEA results (from IterativeEnrichment)
iter_enrich = IterativeEnrichment(...)
igea_results = iter_enrich.results
gene_list_size = len(gene_set.genes)

# Benchmark
benchmark = benchmark_real_results(
    igea_results,
    gene_list_size,
    null_distribution
)

# Interpret results
print(f"Average connections per gene:")
print(f"  Real: {benchmark['real_metrics']['avg_connections_per_gene']:.2f}")
print(f"  Null mean: {benchmark['comparison']['avg_connections_per_gene']['null_mean']:.2f}")
print(f"  Z-score: {benchmark['comparison']['avg_connections_per_gene']['z_score']:.2f}")
print(f"  Percentile: {benchmark['comparison']['avg_connections_per_gene']['percentile']:.1f}%")
```

## Interpretation

### What Good Connectivity Means

**High connectivity** (above null distribution) suggests:
- Genes are functionally related
- Terms are biologically coherent
- The gene list represents a meaningful biological process

**Low connectivity** (at or below null distribution) suggests:
- Genes may be randomly associated
- Terms are not well-connected
- The gene list may not represent a coherent biological process

### Key Metrics to Focus On

1. **avg_connections_per_gene**: 
   - Higher = genes participate in multiple terms (more interconnected)
   - Expected: Real lists > Random lists

2. **network_density**:
   - Higher = more connections relative to possible connections
   - Expected: Real lists > Random lists

3. **n_connected_components**:
   - Lower = more genes/terms are connected (single large component)
   - Expected: Real lists < Random lists (fewer disconnected components)

4. **clustering_coefficient**:
   - Higher = genes in the same terms tend to share other terms
   - Expected: Real lists > Random lists

5. **gene_centrality_max**:
   - Higher = some genes are highly connected (hubs)
   - Expected: Real lists > Random lists

## Example Output

```
Benchmark Results for Gene List (size: 272)

Metric                          Real    Null Mean  Z-score  Percentile
─────────────────────────────────────────────────────────────────────
avg_connections_per_gene        2.45    1.23       3.12     99.9%
network_density                 0.012   0.006      2.89     99.8%
n_connected_components          3       8          -2.45    0.7%
largest_component_size          85      45         4.23     99.99%
gene_centrality_max             12      6          3.67     99.98%
clustering_coefficient          0.34    0.18       2.78     99.7%

Interpretation: Network connectivity is significantly better than random
(p < 0.001 for all metrics). This suggests the gene list represents a
coherent biological process with strong functional relationships.
```

## Advantages of This Approach

1. **Uses existing permutation data**: No need to generate new random lists
2. **Gene list size dependent**: Accounts for the fact that larger lists naturally have more connections
3. **Multiple metrics**: Captures different aspects of connectivity
4. **Statistical rigor**: Provides z-scores and percentiles for interpretation
5. **Hypothesis testing**: Can formally test if connectivity is better than random

## Files

- `code/network_connectivity_benchmark.py`: Main implementation
- `results/connectivity_metrics.tsv`: Connectivity metrics for all permutations
- `results/connectivity_null_distribution.json`: Null distributions by gene list size

## Next Steps

1. Run the benchmark script on permutation data
2. Test on real gene lists (e.g., HIV gene list)
3. Validate that real lists show better connectivity
4. Integrate into iGEA workflow for automatic benchmarking
