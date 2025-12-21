# Density Normalization Proposal

## Problem Statement

Cluster density is not directly comparable across different cluster sizes because:
- **Larger clusters naturally have lower density** due to the constraint that each gene can connect to at most one term per library
- The maximum possible edges formula is: `max_edges = cluster_genes × n_libraries`
- As cluster size increases, the ratio of actual edges to max edges tends to decrease
- Comparing a large real cluster (e.g., 200 nodes, density 0.3) against a small random cluster (e.g., 20 nodes, density 0.4) is unfair

## Current Approach

Currently, we compare density directly:
- Real cluster density vs. null distribution of `largest_cluster_density`
- This compares against the average density of largest clusters from random gene lists
- **Problem**: Random largest clusters are typically much smaller than real largest clusters

## Proposed Solution: Size-Binned Density Comparison

### Approach 1: Size-Binned Null Distribution (Recommended)

**Concept**: Compare density against null distribution for clusters of similar size.

**Implementation**:
1. When building null distribution, bin clusters by size:
   - Small: < 20 nodes
   - Medium: 20-50 nodes
   - Large: 50-100 nodes
   - Very Large: 100-200 nodes
   - Extra Large: > 200 nodes

2. For each size bin, compute density statistics:
   - Mean density
   - Standard deviation
   - Min, max, median

3. When benchmarking:
   - Determine which size bin the real cluster belongs to
   - Compare real cluster density against null distribution for that size bin
   - Compute z-score and percentile within the size bin

**Advantages**:
- Fair comparison: compares like with like
- Statistically sound: accounts for size-dependent density variation
- Maintains interpretability: still shows density, just normalized by size

**Disadvantages**:
- Requires storing additional data (size-binned statistics)
- More complex implementation
- Need to handle edge cases (very small or very large clusters)

### Approach 2: Expected Density Normalization

**Concept**: Compare observed density to expected density for that cluster size.

**Implementation**:
1. Fit a regression model: `density ~ f(cluster_size)` from null distribution
2. Compute expected density for real cluster size: `expected_density = f(real_cluster_size)`
3. Compute normalized density: `normalized_density = observed_density / expected_density`
4. Compare normalized density against null distribution of normalized densities

**Advantages**:
- Single metric (normalized density ratio)
- Accounts for size-dependent trends
- Can use all permutation data (not just largest clusters)

**Disadvantages**:
- Requires fitting a model (may not be linear)
- Less intuitive interpretation
- Model assumptions may not hold

### Approach 3: Percentile Within Size Bins

**Concept**: Compute percentile rank within size bin.

**Implementation**:
1. Bin clusters by size in null distribution
2. For each size bin, compute density percentiles
3. When benchmarking:
   - Find which size bin the real cluster belongs to
   - Compute percentile rank: what % of clusters in that size bin have lower density?

**Advantages**:
- Non-parametric (no distribution assumptions)
- Easy to interpret: "top X% for clusters of this size"
- Robust to outliers

**Disadvantages**:
- Requires sufficient data in each size bin
- Percentile may be less precise than z-score

## Recommended Implementation: Hybrid Approach

Combine **Approach 1 (Size-Binned)** with **Approach 3 (Percentile)**:

1. **Store size-binned statistics** in null distribution:
   ```python
   {
       'gene_list_size': 200,
       'largest_cluster_density': {
           'mean': 0.35,
           'std': 0.15,
           # Size-binned statistics
           'by_size_bin': {
               'small_<20': {'mean': 0.45, 'std': 0.20, 'n': 150},
               'medium_20_50': {'mean': 0.35, 'std': 0.15, 'n': 200},
               'large_50_100': {'mean': 0.25, 'std': 0.12, 'n': 100},
               'very_large_100_200': {'mean': 0.20, 'std': 0.10, 'n': 50},
               'extra_large_>200': {'mean': 0.15, 'std': 0.08, 'n': 10}
           }
       }
   }
   ```

2. **When benchmarking**:
   - Determine size bin for real cluster
   - Compare against size-binned statistics
   - Compute both z-score and percentile within bin

3. **Display both metrics**:
   - Raw density (for reference)
   - Size-normalized z-score (for comparison)
   - Size-normalized percentile (for interpretation)

## Size Bin Definitions

Suggested bins based on typical cluster sizes:
- **Tiny**: < 10 nodes (very small clusters, rare in large gene lists)
- **Small**: 10-20 nodes
- **Medium**: 20-50 nodes (typical random largest clusters)
- **Large**: 50-100 nodes (common in real gene lists)
- **Very Large**: 100-200 nodes (large real clusters)
- **Extra Large**: > 200 nodes (very large real clusters)

Bins can be adjusted based on actual data distribution.

## Implementation Steps

1. **Modify null distribution computation** to include size-binned statistics
2. **Update `benchmark_cluster()`** to use size-binned comparison for density
3. **Add new metric**: `cluster_density_size_normalized_z` and `cluster_density_size_normalized_percentile`
4. **Update UI** to show both raw and normalized density metrics
5. **Update documentation** to explain the normalization

## Example

**Real cluster**:
- Size: 150 nodes
- Density: 0.25

**Null distribution (all sizes)**:
- Mean density: 0.35
- Std: 0.15

**Null distribution (size bin 100-200)**:
- Mean density: 0.20
- Std: 0.10

**Comparison**:
- Raw z-score: (0.25 - 0.35) / 0.15 = -0.67 (worse than average)
- Size-normalized z-score: (0.25 - 0.20) / 0.10 = 0.50 (better than average for this size)
- Size-normalized percentile: ~69% (top 31% for clusters of this size)

**Conclusion**: The real cluster has lower raw density, but **higher density than expected for its size**, indicating good connectivity.

