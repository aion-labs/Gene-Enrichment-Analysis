# Null Distribution Metrics Fix

## Problem Identified

The statistical benchmarking was showing "N/A" for all metrics because:

1. **Metric Name Mismatch**: The metrics list included metric names that were never computed:
   - `'cluster_size'`, `'n_genes'`, `'n_terms'`, `'n_edges'`, `'density'`, `'n_libraries'`
   - These were in the requested metrics list but not being stored in `all_cluster_stats`

2. **Missing Computations**: Several metrics were in the list but not being computed:
   - `'largest_cluster_libraries'` - was in the list but not computed/stored
   - `'fraction_in_largest_cluster'` - was in the list but not computed/stored
   - `'avg_cluster_libraries'` - was in the list but not computed/stored

3. **Result**: When the code tried to compute statistics for these metrics, they weren't in the DataFrame, so they were skipped with warnings, resulting in empty/null statistics.

## Root Cause

In `compute_null_distribution_from_raw_permutations()`:
- The function was only storing a subset of the requested metrics
- The metrics list included names that didn't match what was being computed
- Network-wide metrics (like `fraction_in_largest_cluster`) require calling `analyzer.compute_metrics()` but this wasn't being done

## Fix Applied

1. **Removed incorrect metric names** from the default metrics list:
   - Removed: `'cluster_size'`, `'n_genes'`, `'n_terms'`, `'n_edges'`, `'density'`, `'n_libraries'`
   - These don't make sense as standalone metrics (they should be prefixed with `largest_cluster_*`)

2. **Added missing computations**:
   - Now computes `network_metrics = analyzer.compute_metrics(include_library_diversity=True)` to get network-wide metrics
   - Stores `largest_cluster_libraries` from the largest cluster
   - Stores `avg_cluster_libraries` computed from all clusters
   - Stores `fraction_in_largest_cluster` from network metrics

3. **Handles empty clusters properly**:
   - When there are no clusters, all metrics are set to 0 (including the new ones)
   - This ensures we can still compute a distribution even if some permutations have no clusters

## Expected Behavior After Fix

- All requested metrics should now be computed and stored
- Statistics (mean, std, min, max, median) should be computed for all metrics
- The null distribution should contain all the metrics needed for benchmarking
- No more "N/A" values in the statistical reports

## Testing

To verify the fix works:
1. Re-run the HIV benchmarking: `python benchmark_hiv_connectivity.py`
2. Check that the null distribution contains all expected metrics
3. Verify that the statistical report shows actual values instead of "N/A"

