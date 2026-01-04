# Benchmarking Methodology Confirmation

## Question
**Is the HIV cluster benchmarked against the largest random clusters or all random clusters?**

## Answer: ✅ **Largest Random Clusters Only**

### Process Flow

1. **Filter Permutations** (by p-value and libraries):
   - Filter merged permutation results by p-value threshold (0.01)
   - Filter by selected libraries (Reactome, KEGG, GO BP, GO MF, GO CC)

2. **For Each Permutation**:
   - Build network from filtered enrichment results
   - Identify all clusters (connected components)
   - **Extract ONLY the largest cluster**: `largest_cluster = clusters[0]` (clusters are sorted by size)
   - Store metrics from the largest cluster only

3. **Build Null Distribution**:
   - Collect largest cluster metrics from all permutations
   - Compute statistics (mean, std, min, max, median) for each metric
   - This creates a distribution of "largest cluster" metrics from random gene lists

4. **Benchmark HIV Cluster**:
   - Compare HIV cluster metrics against null distribution
   - All comparisons use `largest_*` metrics from the null distribution

### Code Verification

#### In `compute_null_distribution_from_raw_permutations`:
```python
# Get clusters (sorted by size, largest first)
clusters = analyzer.get_clusters()

# Extract ONLY the largest cluster
largest_cluster = clusters[0]  # Already sorted by size

# Store metrics from largest cluster only
all_cluster_stats.append({
    'largest_cluster_genes': largest_cluster['n_genes'],
    'largest_cluster_terms': largest_cluster['n_terms'],
    'largest_cluster_edges': largest_cluster['n_edges'],
    'largest_cluster_density': largest_cluster['density'],
    'largest_cluster_libraries': largest_cluster.get('n_libraries', 0),
    'largest_component_size': largest_cluster['size'],
    # ... other metrics
})
```

#### In `benchmark_cluster`:
```python
# Map cluster metrics to null distribution metric names
metric_mapping = {
    'cluster_size': 'largest_component_size',      # Largest cluster only
    'cluster_genes': 'largest_cluster_genes',      # Largest cluster only
    'cluster_terms': 'largest_cluster_terms',      # Largest cluster only
    'cluster_edges': 'largest_cluster_edges',     # Largest cluster only
    'cluster_density': 'largest_cluster_density',  # Largest cluster only
    'cluster_libraries': 'largest_cluster_libraries',  # Largest cluster only
}
```

### Confirmation

✅ **The implementation is correct**: 
- After p-value and library filtering
- We extract the largest cluster from each permutation file
- We use only the largest cluster metrics to build the null distribution
- We benchmark each HIV cluster against this null distribution of largest random clusters

### Note on Other Metrics

The code also computes `avg_cluster_*` metrics (using all clusters), but these are **NOT used for benchmarking**. They're computed but not compared. Only `largest_cluster_*` metrics are used in the statistical comparison.

### Result

The HIV cluster benchmarking correctly compares:
- **HIV cluster** (could be any cluster from HIV results)
- **vs. Largest random clusters** (from each permutation)

This is the appropriate comparison because:
- We want to know if the HIV cluster is better than typical largest random clusters
- Comparing against all random clusters (including small ones) would be less meaningful
- The largest cluster is the most informative metric for connectivity

