# Cluster-Based Benchmarking Approach

## Summary

**Key Insight**: In iGEA networks, **clusters (connected components) are the natural unit of connectivity**. Terms can only connect through shared genes, and each term has a fixed gene set from its library. Therefore, we should benchmark at the **cluster level**, not using global network density.

## Why Cluster-Based?

### Problem with Global Network Density

The formula `network_density = n_edges / (n_genes × n_terms)` assumes any gene can connect to any term, but this is **incorrect** for iGEA:

1. **Terms have fixed gene sets**: Each term comes from a library with a predefined set of genes
2. **Not all connections are possible**: A term can only connect to genes that are actually in that term's gene set
3. **Terms connect only through genes**: Terms can only be connected to each other through shared genes

### Correct Approach: Cluster-Level Metrics

**Clusters** (connected components) represent groups of terms that share genes:
- A cluster = a set of terms connected through shared genes
- Cluster density = edges within cluster / (genes in cluster × n_libraries_in_cluster)
- This is correct because in iGEA, each gene can connect to at most 1 term per library

## Benchmarking Strategy

### Tier 1: Primary Cluster Metrics (Focus on Largest Cluster)

These metrics focus on the **largest cluster**, which should be the most informative:

1. **`largest_cluster_genes`** ⭐⭐⭐
   - Number of genes in largest cluster
   - **Hypothesis**: Real lists > Random lists
   - **Interpretation**: More genes participate in the main connected component

2. **`largest_cluster_terms`** ⭐⭐⭐
   - Number of terms in largest cluster
   - **Hypothesis**: Real lists > Random lists
   - **Interpretation**: More terms are connected together

3. **`largest_cluster_edges`** ⭐⭐⭐
   - Number of edges in largest cluster
   - **Hypothesis**: Real lists > Random lists
   - **Interpretation**: More connections within the main cluster

4. **`largest_cluster_density`** ⭐⭐⭐
   - Density = `largest_cluster_edges / (largest_cluster_genes × n_libraries_in_cluster)`
   - In iGEA, each gene can connect to at most 1 term per library
   - **Hypothesis**: Real lists ≥ Random lists
   - **Interpretation**: Higher density = more connections per possible connection

5. **`largest_component_size`** ⭐⭐
   - Total nodes (genes + terms) in largest cluster
   - **Hypothesis**: Real lists > Random lists
   - **Interpretation**: Larger main cluster = better connectivity

6. **`fraction_in_largest_cluster`** ⭐⭐
   - Fraction of all nodes in largest cluster
   - **Hypothesis**: Real lists > Random lists
   - **Interpretation**: More of the network is connected

7. **`fraction_edges_in_largest_cluster`** ⭐⭐
   - Fraction of all edges in largest cluster
   - **Hypothesis**: Real lists > Random lists
   - **Interpretation**: More edges are in the main cluster

### Tier 2: Cluster Distribution Metrics

Characterize the overall cluster structure:

8. **`n_connected_components`** ⭐⭐
   - Number of clusters (lower is better)
   - **Hypothesis**: Real lists < Random lists
   - **Interpretation**: Fewer clusters = more connected network

9. **`avg_cluster_size`** ⭐
   - Average size across all clusters
   - **Hypothesis**: Real lists > Random lists
   - **Interpretation**: Larger average clusters

10. **`weighted_avg_cluster_density`** ⭐
    - Size-weighted average cluster density
    - **Hypothesis**: Real lists ≥ Random lists
    - **Interpretation**: Larger clusters tend to be denser

11. **`avg_cluster_density`** ⭐
    - Average density across all clusters
    - **Hypothesis**: Real lists ≥ Random lists
    - **Interpretation**: All clusters are more connected

### Tier 3: Gene Connectivity Metrics

Measure how well genes connect terms (within clusters):

12. **`hub_genes_count`** ⭐
    - Number of genes connecting to ≥3 terms
    - **Hypothesis**: Real lists > Random lists
    - **Interpretation**: More hub genes = better connectivity

13. **`avg_connections_per_gene`** ⭐
    - Average number of terms each gene connects to
    - **Hypothesis**: Real lists > Random lists
    - **Interpretation**: Genes participate in more terms

14. **`gene_centrality_max`** ⭐
    - Maximum number of terms a single gene connects to
    - **Hypothesis**: Real lists > Random lists
    - **Interpretation**: Presence of highly connected hub genes

### Tier 4: Global Clustering

15. **`clustering_coefficient`** ⭐
    - Bipartite clustering coefficient
    - **Hypothesis**: Real lists ≥ Random lists
    - **Interpretation**: Terms sharing genes tend to share other genes

## Implementation

### Current Status

✅ **Already implemented**:
- Cluster detection (BFS on bipartite graph)
- Cluster metrics (size, density, genes, terms, edges)
- Largest cluster metrics
- Cluster distribution metrics
- Gene connectivity metrics

✅ **Metrics are computed correctly**:
- Cluster density = `cluster_edges / (cluster_genes × n_libraries_in_cluster)`
- In iGEA, each gene can connect to at most 1 term per library
- Only considers actual connections within each cluster

### Benchmarking Process

1. **For each permutation**:
   - Compute all cluster metrics
   - Store per gene list size

2. **Build null distribution**:
   - Stratify by gene list size
   - Compute statistics (mean, std, percentiles) for each metric

3. **For real gene list**:
   - Compute same metrics
   - Compare against null distribution (same gene list size)
   - Calculate z-scores and percentiles
   - Report which metrics are significantly better

### Interpretation

**Real gene lists should show**:
- ✅ **Larger largest cluster** (more genes and terms)
- ✅ **Higher largest cluster density** (more connections)
- ✅ **Fewer clusters** (more connected network)
- ✅ **More hub genes** (genes connecting multiple terms)
- ✅ **Higher fraction in largest cluster** (more of network is connected)

**Example**:
- **Real list**: Largest cluster = 50 genes, 20 terms, 200 edges, density = 0.20
- **Random list**: Largest cluster = 10 genes, 5 terms, 15 edges, density = 0.30
- **Conclusion**: Real list has much larger cluster (better connectivity) even if density is similar

## Summary

**Focus**: Cluster-based metrics, especially the largest cluster
**Key metrics**: `largest_cluster_genes`, `largest_cluster_terms`, `largest_cluster_edges`, `largest_cluster_density`
**Hypothesis**: Real gene lists form fewer, larger, more connected clusters than random lists
