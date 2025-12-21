# Benchmark Workflow Explanation

## Current Problem

**The benchmark script ALWAYS reruns iGEA from scratch** - there is NO caching mechanism. Every execution:
1. Runs full iGEA enrichment for all 5 libraries (takes 2-5 minutes)
2. Computes network metrics
3. Generates reports

This is why commands take so long - Python IS running, but it's doing expensive iGEA computations.

## Current Workflow (benchmark_hiv_connectivity.py)

### Step 1-5: Setup (Fast)
- Load HIV gene list
- Load libraries
- Load background
- Find libraries with permutation data
- Start null distribution computation in parallel thread

### Step 6: Null Distribution (Parallel, Fast)
- Runs in background thread
- Loads Parquet files (pre-computed)
- Filters by libraries
- Computes statistics

### Step 7: iGEA Enrichment (SLOW - 2-5 minutes)
```python
# Lines 369-387: ALWAYS runs iGEA
for lib_name, library in libraries.items():
    iter_enrich = run_igea_for_library(gene_set, library, background)  # SLOW!
    all_results[lib_name] = iter_enrich
```

**This is the bottleneck** - no way to skip it.

### Step 8: Network Analysis (Fast)
- Build network from iGEA results
- Compute clusters
- Benchmark against null distribution
- Generate reports

## Solution: Add Caching

The script should:
1. Check if iGEA results already exist (from previous run)
2. If they exist, load them instead of rerunning
3. Only regenerate reports

## Why Commands Appear to Hang

When you run `generate_cluster_statistical_report.py` directly, it should be fast (just reads TSV and formats). But if you're running `benchmark_hiv_connectivity.py`, it's doing the full iGEA run.

## Quick Fix for Report Regeneration

To regenerate just the report without rerunning iGEA:

```bash
python3 generate_cluster_statistical_report.py \
  results/hiv_connectivity_benchmark/hiv_clusters_by_size.tsv \
  results/hiv_connectivity_benchmark/hiv_clusters_statistical_report.txt \
  "HIV Gene List"
```

This should take < 1 second.

