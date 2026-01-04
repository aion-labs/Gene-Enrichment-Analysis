# Null Distribution Methodology Verification

## Current Implementation Analysis

### Question
**Is the HIV cluster benchmarked against the largest random clusters or all random clusters?**

### Current Behavior

#### For Raw Permutations (`compute_null_distribution_from_raw_permutations`)

1. **After p-value and library filtering**: ✅ Correct
   - Filters permutation data by p-value threshold (0.01)
   - Filters by selected libraries (Reactome, KEGG, GO BP, GO MF, GO CC)

2. **Extract largest cluster from each permutation**: ✅ Partially Correct
   - Gets all clusters: `clusters = analyzer.get_clusters()` (sorted by size, largest first)
   - Extracts largest cluster: `largest_cluster = clusters[0]`
   - Uses largest cluster for:
     - `largest_cluster_genes` ✅
     - `largest_cluster_terms` ✅
     - `largest_cluster_edges` ✅
     - `largest_cluster_density` ✅
     - `largest_cluster_libraries` ✅
     - `largest_component_size` ✅
   
   - BUT also computes metrics using ALL clusters:
     - `avg_cluster_size` ❌ (uses all clusters)
     - `avg_cluster_density` ❌ (uses all clusters)
     - `avg_cluster_libraries` ❌ (uses all clusters)
     - `n_connected_components` ❌ (count of all clusters)

3. **Benchmarking**: ✅ Correct
   - `benchmark_cluster()` only uses metrics that map to `largest_*` metrics
   - The mapping is:
     - `cluster_size` → `largest_component_size`
     - `cluster_genes` → `largest_cluster_genes`
     - `cluster_terms` → `largest_cluster_terms`
     - `cluster_edges` → `largest_cluster_edges`
     - `cluster_density` → `largest_cluster_density`
     - `cluster_libraries` → `largest_cluster_libraries`

#### For Parquet-based Computation (`compute_null_distribution_from_parquet`)

1. **Uses helper function `get_largest_cluster_metric`**: ✅ Correct
   - Specifically extracts `cluster_number == 1` (the largest cluster)
   - All metrics use this function to get only the largest cluster

## Conclusion

**Current Status**: The benchmarking is **correctly** using only the largest cluster from each permutation for the metrics that are actually compared.

However, there's a **minor inconsistency**:
- The code computes `avg_cluster_*` metrics using all clusters, but these are NOT used for benchmarking
- Only `largest_cluster_*` metrics are used for benchmarking

## Recommendation

For clarity and consistency, we should:
1. **Document** that only largest clusters are used for benchmarking
2. **Optionally remove** `avg_cluster_*` metrics from the null distribution computation if they're not used
3. **Or keep them** for potential future use (they might be useful for other analyses)

The current implementation is **functionally correct** for benchmarking purposes, as only the largest cluster metrics are used in the comparison.

