# Cluster-Based Connectivity Benchmarking Proposal

## Problem with Current Approach

The current `network_density` calculation is **incorrect**:
- **Current**: `density = n_edges / (n_genes × n_terms)`
- **Issue**: Assumes any gene can connect to any term
- **Reality**: A term can only connect to genes that are actually in that term's gene set

## Correct Understanding of iGEA Network

1. **Terms come from libraries**: Each term has a predefined gene set
2. **Gene-term edges**: Only exist if the gene is in the term's gene set
3. **Term-term connections**: Terms can only connect through shared genes (not directly)
4. **Clusters**: Connected components where terms share genes

## Proposed Cluster-Based Metrics

### 1. Cluster-Level Metrics (Primary)

For each connected component (cluster):

#### Cluster Size Metrics
- **`cluster_size`**: Total nodes in cluster (genes + terms)
- **`cluster_genes`**: Number of genes in cluster
- **`cluster_terms`**: Number of terms in cluster
- **`cluster_edges`**: Number of edges within cluster

#### Cluster Density (Correct Calculation)
- **`cluster_density`**: For each cluster, calculate:
  ```
  cluster_density = cluster_edges / (cluster_genes × n_libraries_in_cluster)
  ```
  This is correct because in iGEA, each gene can connect to at most 1 term per library.
  So the maximum possible edges = genes × number of libraries in the cluster.

#### Cluster Connectivity
- **`avg_terms_per_gene_in_cluster`**: Average number of terms each gene connects to within the cluster
- **`avg_genes_per_term_in_cluster`**: Average number of genes each term has within the cluster

### 2. Network-Level Summary Metrics

Aggregate across all clusters:

#### Cluster Distribution
- **`n_clusters`**: Total number of connected components
- **`largest_cluster_size`**: Size (genes + terms) of largest cluster
- **`largest_cluster_genes`**: Number of genes in largest cluster
- **`largest_cluster_terms`**: Number of terms in largest cluster
- **`largest_cluster_edges`**: Number of edges in largest cluster
- **`largest_cluster_density`**: Density of largest cluster
- **`avg_cluster_size`**: Average cluster size
- **`median_cluster_size`**: Median cluster size
- **`cluster_size_std`**: Standard deviation of cluster sizes

#### Cluster Quality Metrics
- **`avg_cluster_density`**: Average density across all clusters
- **`weighted_avg_cluster_density`**: Density weighted by cluster size
- **`fraction_in_largest_cluster`**: Fraction of all nodes in largest cluster
- **`fraction_edges_in_largest_cluster`**: Fraction of all edges in largest cluster

### 3. Gene Connectivity Metrics (Secondary)

- **`avg_terms_per_gene`**: Average number of terms each gene connects to (already computed)
- **`gene_centrality_distribution`**: Distribution of gene degrees
- **`hub_genes_count`**: Number of genes connecting to ≥3 terms

### 4. Term Connectivity Metrics (Secondary)

- **`avg_genes_per_term`**: Average number of genes in each term (already computed)
- **`term_centrality_distribution`**: Distribution of term sizes

## Benchmarking Strategy

### Per-Cluster Benchmarking

For each cluster size category (e.g., small: 2-5 nodes, medium: 6-20, large: 21+):

1. **Extract clusters** from permutation results
2. **Stratify by cluster size** (and gene list size)
3. **Build null distributions** for:
   - Cluster density
   - Cluster size distribution
   - Number of clusters per network
   - Largest cluster metrics

### Network-Level Benchmarking

Compare aggregate metrics:
- Number of clusters (fewer = better connectivity)
- Largest cluster size (larger = better connectivity)
- Largest cluster density (higher = better connectivity)
- Fraction of nodes in largest cluster (higher = better connectivity)

## Implementation Plan

### Step 1: Enhance Cluster Detection

Modify `_compute_connected_components()` to return:
- Full cluster information (genes, terms, edges)
- Cluster-level metrics

### Step 2: Add Cluster Metrics

Add methods to compute:
- Cluster density (correct calculation)
- Cluster size distribution
- Largest cluster metrics

### Step 3: Update Benchmarking

- Remove incorrect `network_density`
- Add cluster-based metrics to null distribution
- Benchmark by cluster size categories

### Step 4: Update Documentation

- Explain cluster-based approach
- Document correct density calculation
- Update interpretation guidelines

## Expected Results

For real gene lists vs. random:
- **Fewer clusters**: More terms connected through shared genes
- **Larger largest cluster**: More genes and terms in main component
- **Higher cluster density**: More connections within clusters
- **Higher fraction in largest cluster**: More network cohesion

## Example Interpretation

**Good connectivity** (real gene list):
- 1-2 large clusters containing most genes/terms
- High cluster density (>0.1)
- Large fraction of nodes in largest cluster (>0.7)

**Poor connectivity** (random gene list):
- Many small clusters (5-10+)
- Low cluster density (<0.05)
- Small fraction of nodes in largest cluster (<0.3)
