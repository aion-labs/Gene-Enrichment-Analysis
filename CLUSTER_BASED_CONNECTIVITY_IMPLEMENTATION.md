# Cluster-Based Connectivity Implementation

## Changes Made

### 1. Removed Incorrect Network Density

**Problem**: The original `network_density = n_edges / (n_genes × n_terms)` assumed any gene could connect to any term, which is incorrect.

**Solution**: Removed this metric entirely. Density is now computed **per cluster** where it's meaningful.

### 2. Enhanced Cluster Detection

**Before**: `_compute_connected_components()` only returned count and largest size.

**After**: Now returns full cluster information:
- Genes and terms in each cluster
- Number of edges within each cluster
- Cluster density (correct calculation: `edges / (cluster_genes × n_libraries_in_cluster)`)
- In iGEA, each gene can connect to at most 1 term per library

### 3. New Cluster-Based Metrics

#### Primary Metrics (for benchmarking)

**Cluster Size Metrics**:
- `largest_cluster_genes`: Number of genes in largest cluster
- `largest_cluster_terms`: Number of terms in largest cluster
- `largest_cluster_edges`: Number of edges in largest cluster
- `largest_cluster_density`: Density of largest cluster (correct calculation)
- `avg_cluster_size`: Average cluster size
- `median_cluster_size`: Median cluster size
- `cluster_size_std`: Standard deviation of cluster sizes

**Cluster Density Metrics**:
- `avg_cluster_density`: Average density across all clusters
- `weighted_avg_cluster_density`: Density weighted by cluster size

**Cluster Connectivity Metrics**:
- `fraction_in_largest_cluster`: Fraction of all nodes in largest cluster
- `fraction_edges_in_largest_cluster`: Fraction of all edges in largest cluster
- `n_connected_components`: Number of clusters (lower is better)

#### Secondary Metrics

**Gene Connectivity**:
- `avg_connections_per_gene`: Average terms per gene
- `gene_centrality_max`: Maximum terms any gene connects to
- `hub_genes_count`: Number of genes connecting to ≥3 terms

**Global Clustering**:
- `clustering_coefficient`: Bipartite clustering coefficient

### 4. Updated Benchmarking Logic

**"Is Better" Logic**:
- `n_connected_components`: **Lower is better** (fewer clusters = better connectivity)
- All other metrics: **Higher is better**

## Interpretation

### Good Connectivity (Real Gene Lists)

- **Few clusters**: `n_connected_components` = 1-2
- **Large largest cluster**: `largest_component_size` > 50% of network
- **High cluster density**: `largest_cluster_density` > 0.1
- **High fraction in largest**: `fraction_in_largest_cluster` > 0.7
- **Many hub genes**: `hub_genes_count` > 5

### Poor Connectivity (Random Gene Lists)

- **Many small clusters**: `n_connected_components` = 5-10+
- **Small largest cluster**: `largest_component_size` < 30% of network
- **Low cluster density**: `largest_cluster_density` < 0.05
- **Low fraction in largest**: `fraction_in_largest_cluster` < 0.3
- **Few hub genes**: `hub_genes_count` < 2

## Next Steps

1. **Rebuild null distribution** with new cluster-based metrics
2. **Re-run HIV benchmarking** to see cluster-based results
3. **Compare** cluster metrics vs. old metrics

## Example Usage

```python
from code.network_connectivity_benchmark import NetworkConnectivityAnalyzer

analyzer = NetworkConnectivityAnalyzer()
analyzer.add_igea_results(igea_results)
metrics = analyzer.compute_metrics()

# Key cluster metrics
print(f"Number of clusters: {metrics['n_connected_components']}")
print(f"Largest cluster: {metrics['largest_component_size']} nodes")
print(f"Largest cluster density: {metrics['largest_cluster_density']:.3f}")
print(f"Fraction in largest: {metrics['fraction_in_largest_cluster']:.2%}")
```
