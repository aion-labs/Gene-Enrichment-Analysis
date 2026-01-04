# Network Connectivity Benchmarking - Summary

## Strategy Overview

**Hypothesis**: Real gene lists will have **better network connectivity** than random gene lists, as measured by graph metrics on the iGEA network.

**Approach**: Use permutation results (random gene lists) to establish a **null distribution** of connectivity metrics, stratified by gene list size.

## Network Structure

iGEA creates a **bipartite graph**:
- **Nodes**: Genes (from input list) + Terms (enriched pathways)
- **Edges**: Gene-term connections (gene is a member of that term)

## Connectivity Metrics

### Primary Metrics

1. **avg_connections_per_gene** ⭐
   - Average number of terms each gene connects to
   - **Expected**: Real lists > Random lists
   - **Interpretation**: Higher = genes participate in multiple pathways (more interconnected)

2. **network_density** ⭐
   - Ratio of actual edges to possible edges
   - **Expected**: Real lists > Random lists
   - **Interpretation**: Higher = more connections relative to network size

3. **n_connected_components**
   - Number of disconnected subgraphs
   - **Expected**: Real lists < Random lists
   - **Interpretation**: Lower = more genes/terms are connected (single large component)

4. **clustering_coefficient** ⭐
   - Measure of how interconnected the network is
   - **Expected**: Real lists > Random lists
   - **Interpretation**: Higher = genes in same terms tend to share other terms

5. **gene_centrality_max**
   - Maximum number of terms a single gene connects to
   - **Expected**: Real lists > Random lists
   - **Interpretation**: Higher = presence of highly connected hub genes

### Secondary Metrics

- **largest_component_size**: Size of largest connected component
- **gene_centrality_mean/std**: Distribution of gene degrees
- **term_centrality_mean/std/max**: Distribution of term sizes

## Implementation

### Files Created

1. **`code/network_connectivity_benchmark.py`**
   - `NetworkConnectivityAnalyzer`: Computes connectivity metrics from iGEA results
   - `compute_connectivity_from_permutation_results()`: Processes permutation data
   - `build_null_distribution()`: Creates null distributions by gene list size
   - `benchmark_real_results()`: Compares real results against null

2. **`example_network_benchmark.py`**
   - Example usage script
   - Demonstrates computing null distribution and benchmarking

3. **`NETWORK_CONNECTIVITY_BENCHMARKING.md`**
   - Detailed documentation

## Usage Workflow

### Step 1: Compute Null Distribution (One-Time)

```bash
python code/network_connectivity_benchmark.py \
    --permutation-file archive/permutation_analysis/results/permutation_results/merged_permutation_results.tsv \
    --output-metrics results/connectivity_metrics.tsv \
    --output-null results/connectivity_null_distribution.json
```

**Outputs**:
- `connectivity_metrics.tsv`: Metrics for all 20,000 permutations
- `connectivity_null_distribution.json`: Null distributions by gene list size

### Step 2: Benchmark Real Results

```python
from code.network_connectivity_benchmark import benchmark_real_results
import json

# Load null distribution
with open('results/connectivity_null_distribution.json') as f:
    null_dist = json.load(f)

# Get iGEA results
iter_enrich = IterativeEnrichment(...)
results = iter_enrich.results
gene_list_size = len(gene_set.genes)

# Benchmark
benchmark = benchmark_real_results(results, gene_list_size, null_dist)

# Interpret
for metric, comp in benchmark['comparison'].items():
    print(f"{metric}: z-score={comp['z_score']:.2f}, percentile={comp['percentile']:.1f}%")
```

## Advantages

✅ **Uses existing permutation data** - No need to generate new random lists  
✅ **Gene list size dependent** - Accounts for natural scaling  
✅ **Multiple metrics** - Captures different aspects of connectivity  
✅ **Statistical rigor** - Z-scores and percentiles for interpretation  
✅ **Hypothesis testing** - Can formally test if connectivity > random  

## Expected Results

For **real gene lists** (e.g., HIV, disease signatures):
- **avg_connections_per_gene**: Z-score > 2 (top 2.5%)
- **network_density**: Z-score > 2
- **clustering_coefficient**: Z-score > 2
- **n_connected_components**: Z-score < -2 (fewer components = better)

For **random gene lists**:
- All metrics should be near null mean (Z-score ≈ 0)
- Percentiles around 50%

## Next Steps

1. ✅ **Implementation complete** - Scripts are ready
2. ⏳ **Run on permutation data** - Compute null distributions (takes ~5-10 minutes)
3. ⏳ **Test on real lists** - Validate hypothesis (e.g., HIV gene list)
4. ⏳ **Integrate into iGEA** - Add automatic benchmarking to workflow

## Files

- `code/network_connectivity_benchmark.py` - Main implementation
- `example_network_benchmark.py` - Usage examples
- `NETWORK_CONNECTIVITY_BENCHMARKING.md` - Detailed documentation
- `NETWORK_BENCHMARK_SUMMARY.md` - This summary
