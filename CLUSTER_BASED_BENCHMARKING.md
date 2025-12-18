# Cluster-Based Connectivity Benchmarking

## Problem with Global Network Density

The previous `network_density = n_edges / (n_genes × n_terms)` calculation was **incorrect** because:

1. **Terms have fixed gene sets**: Each term comes from a library and can only connect to genes that are actually in that term's predefined gene set
2. **Not all gene-term pairs are possible**: A term cannot connect to any arbitrary gene
3. **Terms connect only through genes**: Terms can only be connected to each other through shared genes, not directly

## Correct Understanding of iGEA Network

### Network Structure
- **Bipartite graph**: Genes and terms are separate node types
- **Edges**: Only exist if a gene is actually in a term's gene set (from the library)
- **Clusters**: Connected components where terms share genes
- **Connectivity**: Measured by how well genes connect multiple terms, and how terms cluster together

### Key Insight
**Clusters are the natural unit of connectivity** in iGEA:
- A cluster represents a group of terms that share genes
- Larger, denser clusters indicate better biological coherence
- Real gene lists should form fewer, larger, more connected clusters than random lists

## Proposed Cluster-Based Benchmarking Approach

### 1. Primary Metrics: Largest Cluster

Focus on the **largest cluster** as the primary indicator of connectivity:

#### Size Metrics
- **`largest_cluster_genes`**: Number of genes in largest cluster
- **`largest_cluster_terms`**: Number of terms in largest cluster  
- **`largest_cluster_size`**: Total nodes (genes + terms) in largest cluster
- **`largest_cluster_edges`**: Number of edges within largest cluster

#### Density Metrics (Correct Calculation)
- **`largest_cluster_density`**: 
  ```
  density = largest_cluster_edges / (largest_cluster_genes × n_libraries_in_cluster)
  ```
  This is correct because in iGEA, each gene can connect to at most 1 term per library.
  So the maximum possible edges = genes × number of libraries in the cluster.

#### Coverage Metrics
- **`fraction_in_largest_cluster`**: What fraction of all nodes are in the largest cluster
- **`fraction_edges_in_largest_cluster`**: What fraction of all edges are in the largest cluster

**Hypothesis**: Real gene lists will have:
- Larger largest cluster (more genes and terms)
- Higher largest cluster density (more connections per possible connection)
- Higher fraction in largest cluster (more of the network is connected)

### 2. Cluster Distribution Metrics

Characterize the overall cluster structure:

#### Cluster Count
- **`n_connected_components`**: Total number of clusters
- **Expected**: Real lists < Random lists (fewer, larger clusters)

#### Cluster Size Distribution
- **`avg_cluster_size`**: Average size across all clusters
- **`median_cluster_size`**: Median cluster size
- **`cluster_size_std`**: Standard deviation of cluster sizes
- **Expected**: Real lists have larger average/median cluster sizes

#### Cluster Density Distribution
- **`avg_cluster_density`**: Average density across all clusters
- **`weighted_avg_cluster_density`**: Average density weighted by cluster size
- **Expected**: Real lists have higher cluster densities

### 3. Gene-Level Connectivity (Within Clusters)

Measure how well genes connect terms:

- **`avg_connections_per_gene`**: Average number of terms each gene connects to
- **`gene_centrality_max`**: Maximum number of terms a single gene connects to
- **`hub_genes_count`**: Number of genes connecting to ≥3 terms
- **Expected**: Real lists have more hub genes and higher average connections

### 4. Benchmarking Strategy

#### Option A: Aggregate Cluster Metrics (Current Approach)
- Compute cluster metrics for each permutation
- Build null distribution of aggregate metrics (e.g., largest_cluster_size, largest_cluster_density)
- Compare real results against null distribution

**Pros**: Simple, single comparison
**Cons**: Loses information about cluster distribution

#### Option B: Cluster-by-Cluster Comparison (Recommended)
For each permutation and real result:
1. Identify all clusters
2. Sort clusters by size (largest first)
3. Compare top N clusters against null distribution

**Metrics to benchmark**:
- **Largest cluster**: size, density, genes, terms, edges
- **2nd largest cluster**: size, density (if exists)
- **3rd largest cluster**: size, density (if exists)
- **Cluster distribution**: number of clusters, size distribution

**Pros**: More detailed, captures cluster structure
**Cons**: More complex, need to handle variable number of clusters

#### Option C: Hybrid Approach (Recommended)
Combine both:
1. **Primary**: Aggregate metrics (largest cluster size, density, etc.)
2. **Secondary**: Cluster distribution (number of clusters, size distribution)
3. **Tertiary**: Gene-level connectivity (hub genes, average connections)

## Implementation Recommendations

### Metrics to Benchmark (Priority Order)

#### Tier 1: Primary Cluster Metrics
1. **`largest_cluster_size`** ⭐⭐⭐
   - Total nodes in largest cluster
   - **Expected**: Real >> Random
   
2. **`largest_cluster_density`** ⭐⭐⭐
   - Density of largest cluster
   - **Expected**: Real > Random
   
3. **`largest_cluster_terms`** ⭐⭐
   - Number of terms in largest cluster
   - **Expected**: Real > Random

4. **`fraction_in_largest_cluster`** ⭐⭐
   - What fraction of network is in largest cluster
   - **Expected**: Real > Random

#### Tier 2: Cluster Distribution
5. **`n_connected_components`** ⭐⭐
   - Number of clusters (lower is better)
   - **Expected**: Real < Random

6. **`avg_cluster_size`** ⭐
   - Average cluster size
   - **Expected**: Real > Random

7. **`weighted_avg_cluster_density`** ⭐
   - Size-weighted average cluster density
   - **Expected**: Real > Random

#### Tier 3: Gene Connectivity
8. **`hub_genes_count`** ⭐
   - Genes connecting to ≥3 terms
   - **Expected**: Real > Random

9. **`avg_connections_per_gene`** ⭐
   - Average terms per gene
   - **Expected**: Real > Random

### Null Distribution Building

For each permutation:
1. Compute all cluster metrics
2. Extract largest cluster metrics
3. Compute cluster distribution metrics
4. Store per gene list size

### Benchmarking Real Results

1. Compute same metrics for real gene list
2. Compare against null distribution (same gene list size)
3. Calculate z-scores and percentiles
4. Report which metrics are significantly better/worse

## Example Interpretation

**Real gene list (HIV)**:
- Largest cluster: 50 genes, 20 terms, 200 edges
- Largest cluster density: 0.20 (200 / (50 × 20))
- 2 clusters total
- 80% of network in largest cluster

**Random gene list (null)**:
- Largest cluster: 10 genes, 5 terms, 15 edges
- Largest cluster density: 0.30 (15 / (10 × 5)) - but much smaller!
- 15 clusters total
- 20% of network in largest cluster

**Conclusion**: Real list has much larger clusters (better connectivity) even if individual cluster density might be similar.

## Summary

**Focus on clusters, not global network density**:
- ✅ Cluster size (genes + terms in largest cluster)
- ✅ Cluster density (within-cluster connectivity)
- ✅ Cluster distribution (number and sizes)
- ✅ Gene connectivity (how many terms per gene)
- ❌ Global network density (incorrect for this use case)
