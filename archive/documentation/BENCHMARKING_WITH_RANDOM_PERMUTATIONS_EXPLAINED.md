# Benchmarking Process Using Random Permutations - Explained

## Overview

The benchmarking process in iGEA uses **random permutations** to establish a statistical null distribution that allows us to determine whether a real gene list's network connectivity is significantly better than what would be expected by chance. This is a form of **permutation testing** - a non-parametric statistical method.

## The Core Hypothesis

**Real biological gene lists should form better-connected networks than random gene lists of the same size.**

This hypothesis is based on the biological principle that genes involved in related processes (e.g., a disease pathway) work together and share functional relationships, which should be reflected in their enrichment patterns and network connectivity.

---

## The Complete Process

### Phase 1: Generating Random Permutations (Pre-computation)

This phase is done **once** to create a reference dataset. It's computationally expensive but only needs to be done once.

#### Step 1.1: Define the Background
- Start with all known genes in the organism (e.g., ~20,000 human genes)
- This serves as the "pool" from which random genes are sampled
- The background is stored in `data/backgrounds/all_genes.txt`

#### Step 1.2: Generate Random Gene Lists
For each permutation:
1. **Randomly sample** a specific number of genes from the background
   - Sizes tested: 50, 100, 150, 200, 250, 300, 400, 500, 750, 1000 genes
   - Each size gets 1,000 permutations (total ~20,000 permutations)
   - Uses a seed based on permutation index for reproducibility

2. **No biological meaning** - pure random selection
   - Example: `random_200_genes_001`, `random_200_genes_002`, etc.

#### Step 1.3: Run iGEA on Each Random List
Each random gene list is processed through iGEA **exactly like a real gene list**:
- Uses the same libraries (Reactome, KEGG, GO BP, GO MF, GO CC)
- Uses the same parameters (p-value threshold 0.05)
- Runs iterative enrichment analysis
- Generates enrichment results (which pathways are enriched, which genes are in each pathway)

**Implementation**: `code/generate_permutation_distribution.py`
- Runs in parallel (multiprocessing) for efficiency
- Saves results to `results/permutation_results/size_{N}/permutation_{idx}.tsv`
- Can resume if interrupted

#### Step 1.4: Store Permutation Results
Results are stored in two formats:

1. **Individual TSV files** (one per permutation)
   - Location: `results/permutation_results/size_{N}/permutation_{idx}.tsv`
   - Contains: Term, Library, Iteration, p-value, Genes removed, etc.

2. **Merged TSV file** (all permutations combined)
   - Location: `results/permutation_results/merged_permutation_results.tsv`
   - Contains all permutation results in one file
   - Used for flexible p-value filtering

---

### Phase 2: Computing Network Metrics from Permutations

This phase analyzes the network structure of each permutation's enrichment results.

#### Step 2.1: Build Network from iGEA Results
For each permutation, the enrichment results form a **bipartite graph**:
- **Gene nodes**: Individual genes from the random gene list
- **Term nodes**: Enriched biological terms/pathways
- **Edges**: Gene-term connections (a gene is a member of that term)

**Visual Example**:
```
Random Gene List Network:
    GeneX ── Pathway D
    GeneY ── Pathway E  
    GeneZ ── Pathway F
    
This forms THREE SMALL CLUSTERS (no connections between pathways).
```

#### Step 2.2: Identify Clusters (Connected Components)
- **Clusters** are groups of pathways connected through shared genes
- If Pathway A and Pathway B both contain Gene X, they're in the same cluster
- Uses BFS (Breadth-First Search) to find connected components

#### Step 2.3: Compute Connectivity Metrics
For each permutation, compute metrics that measure network connectivity:

**Cluster-Level Metrics** (most important):
- `largest_cluster_genes`: Number of genes in the largest cluster
- `largest_cluster_terms`: Number of pathways in the largest cluster
- `largest_cluster_edges`: Number of gene-pathway connections in largest cluster
- `largest_cluster_density`: Density of largest cluster (connections / max possible)
- `largest_cluster_libraries`: Number of different libraries represented in largest cluster

**Network-Wide Metrics**:
- `n_connected_components`: Number of separate clusters (fewer = better)
- `avg_cluster_size`: Average size across all clusters
- `avg_cluster_density`: Average density across all clusters
- `fraction_in_largest_cluster`: What percentage of genes+terms are in the largest cluster

**Implementation**: `code/network_connectivity_benchmark.py`
- `NetworkConnectivityAnalyzer` class builds the network and computes metrics
- Processes each permutation file and extracts metrics

---

### Phase 3: Building the Null Distribution

This phase aggregates all permutation metrics to create statistical reference distributions.

#### Step 3.1: Stratify by Gene List Size
- Group permutations by gene list size (50, 100, 150, ..., 1000)
- Each size has its own null distribution (because connectivity depends on size)

#### Step 3.2: Compute Statistics for Each Metric
For each gene list size and each metric, compute:
- **Mean** (average value)
- **Standard deviation** (spread of values)
- **Min/Max** (range)
- **Median** (middle value)
- **Percentiles** (25th, 50th, 75th, 90th, 95th, 99th)

**Example Null Distribution for Size 200**:
```
Metric: largest_cluster_genes
  Mean: 64.7 genes
  Std Dev: 20.8
  Min: 10
  Max: 114
  Median: 67
  Percentiles: [45, 67, 82, 95, 102, 108]
```

#### Step 3.3: Store in Efficient Format
Results are stored in **Parquet files** (columnar format, fast to query):
- Location: `results/permutation_cluster_statistics_parquet/cluster_stats_size_{N}.parquet`
- One file per gene list size
- Contains cluster statistics for all permutations of that size
- Organized by library combinations (which libraries were used)

**Implementation**: 
- `code/network_connectivity_benchmark.py` - computes metrics from permutations
- Parquet files are pre-computed and stored for fast access

---

### Phase 4: Runtime Benchmarking (When Analyzing Real Gene Lists)

This phase happens **during** the analysis of a real gene list.

#### Step 4.1: Run iGEA on Real Gene List (Main Thread)
- User submits a real gene list (e.g., HIV-related genes)
- Run iGEA enrichment analysis (same as permutations)
- Get enrichment results

#### Step 4.2: Load Null Distribution (Parallel Thread)
**While enrichment is running**, a separate thread:
1. **Checks which libraries have permutation data**
   - Not all libraries may have been included in permutation generation
   - Only libraries with data can be used for statistics

2. **Loads the appropriate null distribution**
   - Finds Parquet file for the gene list size (or nearest size)
   - Filters to only include libraries that:
     - Were selected by the user
     - Have permutation data available
   
3. **Handles p-value filtering** (if user's threshold < 0.05):
   - If user uses p-value threshold < 0.05, need to filter raw permutation data
   - Loads merged permutation TSV file
   - Filters by user's p-value threshold
   - Recomputes cluster statistics on filtered data
   - This is slower but more accurate

4. **Computes null distribution statistics**
   - Aggregates cluster statistics from Parquet
   - Computes mean, std dev, percentiles for each metric

**Implementation**: `code/parallel_null_distribution.py`
- `compute_null_distribution_from_parquet()` - fast path (uses pre-computed Parquet)
- `compute_null_distribution_from_raw_permutations()` - slow path (filters by p-value)

#### Step 4.3: Build Network from Real Results
After enrichment completes:
1. Build the network from real iGEA results (same as permutations)
2. Find clusters (connected components)
3. Compute the same metrics as for permutations

#### Step 4.4: Statistical Comparison
Compare real metrics to the null distribution:

For each metric:
1. **Get real value** (e.g., `largest_cluster_genes = 174`)
2. **Get null statistics** (e.g., mean = 64.7, std = 20.8)
3. **Compute Z-score**:
   ```
   Z = (real_value - null_mean) / null_std
   Z = (174 - 64.7) / 20.8 = 5.25
   ```
4. **Compute percentile**: What percentage of random lists had values lower than yours?
   - Z = 5.25 → percentile ≈ 100% (top 0.0001%)
5. **Determine significance**:
   - Z > 2.0 and Percentile > 95% → ✓✓ SIGNIFICANTLY BETTER
   - Z > 1.0 and Percentile > 84% → ✓ BETTER
   - -1.0 < Z < 1.0 → ~ SIMILAR
   - Z < -1.0 → ✗ WORSE

**Implementation**: `code/network_connectivity_benchmark.py`
- `benchmark_cluster()` - compares a cluster to null distribution
- `benchmark_real_results()` - benchmarks full network

---

## Key Implementation Details

### File Structure
```
results/
├── permutation_results/
│   ├── size_50/
│   │   ├── permutation_001.tsv
│   │   ├── permutation_002.tsv
│   │   └── ...
│   ├── size_100/
│   │   └── ...
│   └── merged_permutation_results.tsv
│
└── permutation_cluster_statistics_parquet/
    ├── cluster_stats_size_50.parquet
    ├── cluster_stats_size_100.parquet
    └── ...
```

### Size Matching and Interpolation
- If real gene list size doesn't match exactly (e.g., 175 genes):
  - **Round up** to next increment of 50 (175 → 200)
  - Or use **interpolation** between two sizes (e.g., 175 → interpolate between 150 and 200)
  - If size > 1000, cap at 1000 (max available)

### Library Filtering
- Only libraries with permutation data can be used for statistics
- User may select libraries without permutation data → they're included in enrichment but excluded from statistics
- Statistics network uses only libraries with permutation data
- Full network results include all selected libraries

### P-value Threshold Handling
- Permutations were generated with p-value threshold **0.05** (captures all enrichments)
- If user uses p-value threshold ≤ 0.05:
  - Can filter raw permutation data to match user's threshold
  - More accurate but slower (requires reprocessing permutations)
- If user uses p-value threshold > 0.05:
  - Statistical benchmarking is **not available**
  - User gets enrichment results but no statistical comparison

---

## Example Workflow

### Generating Permutations (One-time setup)
```bash
python code/generate_permutation_distribution.py \
    --n-permutations 1000 \
    --n-jobs 32
```

### Benchmarking a Real Gene List
```python
# In benchmark_hiv_connectivity.py or similar

# 1. Load real gene list
gene_set = GeneSet.from_file("hiv_genes.txt")

# 2. Run iGEA (main thread)
iter_enrich = run_igea_for_library(gene_set, library, background)

# 3. Load null distribution (parallel thread)
null_distribution = compute_null_distribution_from_parquet(
    parquet_dir=PARQUET_DIR,
    gene_list_size=gene_set.size,
    selected_libraries=["Reactome", "KEGG", "GO BP"],
    user_p_threshold=0.05
)

# 4. Build network from real results
analyzer = NetworkConnectivityAnalyzer()
analyzer.add_igea_results(iter_enrich.results)
clusters = analyzer.get_clusters()

# 5. Benchmark
benchmark = benchmark_cluster(
    cluster=clusters[0],  # Largest cluster
    null_distribution=null_distribution,
    gene_list_size=gene_set.size
)

# 6. Interpret
print(f"Z-score: {benchmark['z_score']:.2f}")
print(f"Percentile: {benchmark['percentile']:.1f}%")
print(f"Status: {benchmark['status']}")
```

---

## Why This Approach Works

1. **Non-parametric**: Doesn't assume any distribution shape
2. **Gene list size dependent**: Accounts for the fact that larger lists naturally have more connections
3. **Multiple metrics**: Captures different aspects of connectivity (size, density, clustering)
4. **Statistically rigorous**: Provides z-scores and percentiles for interpretation
5. **Fast at runtime**: Pre-computed null distributions load quickly (1-2 seconds)
6. **Flexible**: Can filter by p-value threshold, library selection, etc.

---

## Interpretation Guide

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

1. **largest_cluster_genes**: Higher = more genes working together
2. **largest_cluster_density**: Higher = tighter connections
3. **n_connected_components**: Lower = fewer disconnected groups (better)
4. **fraction_in_largest_cluster**: Higher = more of the network is connected

---

## Summary

The benchmarking process uses random permutations to establish what "typical" network connectivity looks like for random gene lists. By comparing a real gene list's connectivity to this null distribution, we can statistically determine whether the real list shows evidence of biological coherence beyond what would be expected by chance.

The process is:
1. **Pre-compute**: Generate thousands of random gene lists, run iGEA, compute metrics
2. **Aggregate**: Build null distributions stratified by gene list size
3. **Runtime**: Compare real gene list to null distribution (fast, parallel)
4. **Interpret**: Use z-scores and percentiles to assess significance

This provides a rigorous, quantitative way to assess whether gene enrichment results are biologically meaningful or could have occurred by chance.

