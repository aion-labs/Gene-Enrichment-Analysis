# How Density is Computed in iGEA Network Analysis

## Overview

Density is computed **per cluster** (connected component), not at the network level. This is because in iGEA networks, not all gene-term connections are possible - each term has a fixed gene set from its library.

## Formula

For each cluster:

```
cluster_density = cluster_edges / max_possible_edges
```

Where:
- **`cluster_edges`**: Actual number of gene-term connections within the cluster
- **`max_possible_edges`**: Maximum possible edges in the cluster

## Maximum Possible Edges Calculation

The key insight is that **in iGEA, each gene can connect to at most 1 term per library**.

Therefore:
```
max_possible_edges = cluster_genes × n_libraries_in_cluster
```

Where:
- **`cluster_genes`**: Number of unique genes in the cluster
- **`n_libraries_in_cluster`**: Number of different libraries represented in the cluster

## Why This Formula?

### The Constraint

In iGEA's iterative enrichment process:
1. Each term comes from a specific library (e.g., Reactome, KEGG, GO BP)
2. Each gene can only be "removed" (assigned to) **one term per library** per iteration
3. Therefore, a gene can connect to at most **one term per library**

### Example

If a cluster has:
- 10 genes
- 3 libraries (Reactome, KEGG, GO BP)
- Maximum possible edges = 10 × 3 = 30

If the cluster actually has 15 edges, then:
- Density = 15 / 30 = 0.5 (50% of maximum possible connections)

## Implementation Details

The density is computed in `_compute_connected_components()` method:

```python
# Count edges within this cluster
cluster_edges = sum(
    len([t for t in self.gene_to_terms[g] if t in cluster_terms])
    for g in cluster_genes
)

# Count unique libraries in this cluster
cluster_libraries = set()
for term in cluster_terms:
    if term in self.term_to_library:
        cluster_libraries.add(self.term_to_library[term])
n_libraries = len(cluster_libraries)

# Density calculation
max_possible_edges = len(cluster_genes) * n_libraries if n_libraries > 0 else 0
cluster_density = cluster_edges / max_possible_edges if max_possible_edges > 0 else 0.0
```

## Density Metrics Used

1. **`largest_cluster_density`**: Density of the largest cluster
2. **`avg_cluster_density`**: Average density across all clusters
3. **`weighted_avg_cluster_density`**: Density weighted by cluster size

## Interpretation

- **Density = 1.0**: Every gene is connected to a term from every library in the cluster (maximum connectivity)
- **Density = 0.5**: Half of the maximum possible connections exist
- **Density = 0.0**: No connections (shouldn't happen if it's a valid cluster)

## Why Not Global Network Density?

The old formula `network_density = n_edges / (n_genes × n_terms)` was **incorrect** because:
- It assumes any gene can connect to any term
- But terms have fixed gene sets from libraries
- Not all gene-term pairs are possible

Cluster-level density is correct because:
- It accounts for the library constraint
- It measures connectivity within biologically meaningful groups (clusters)
- It reflects how well genes connect terms within a functional group

