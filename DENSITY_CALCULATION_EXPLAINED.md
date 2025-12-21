# How Density is Computed in iGEA Network Analysis

## Overview

Density is computed **per cluster** (connected component) as **average gene connectivity**. This is a simple and intuitive measure: the average number of connections per gene.

## Formula

For each cluster:

```
cluster_density = cluster_edges / n_genes
```

Where:
- **`cluster_edges`**: Actual number of gene-term connections within the cluster
- **`n_genes`**: Number of unique genes in the cluster

## Why This Formula?

### Simple and Intuitive

This formula represents the **average number of connections per gene** in the cluster. It's:
- **Simple**: Easy to understand and interpret
- **Intuitive**: Directly measures how well-connected genes are
- **Comparable**: Can be compared across clusters (though larger clusters may naturally have different connectivity patterns)

### Example

If a cluster has:
- 10 genes
- 15 edges total

Then:
- Density = 15 / 10 = 1.5

This means: On average, each gene connects to 1.5 terms.

## Why This Formula?

### The Constraint

In iGEA's iterative enrichment process:
1. Each term comes from a specific library (e.g., Reactome, KEGG, GO BP)
2. Each gene can only be "removed" (assigned to) **one term per library** per iteration
3. Therefore, a gene can connect to at most **one term per library**

### Example

If a cluster has:
- 10 genes
- 15 edges total

Then:
- Density = 15 / 10 = 1.5

This means: On average, each gene connects to 1.5 terms.

## Implementation Details

The density is computed in `_compute_connected_components()` method:

```python
# Count edges within this cluster
cluster_edges = sum(
    len([t for t in self.gene_to_terms[g] if t in cluster_terms])
    for g in cluster_genes
)

# Density calculation: Average gene connectivity
if len(cluster_genes) > 0:
    cluster_density = cluster_edges / len(cluster_genes)
else:
    cluster_density = 0.0
```

## Density Metrics Used

1. **`largest_cluster_density`**: Density of the largest cluster
2. **`avg_cluster_density`**: Average density across all clusters
3. **`weighted_avg_cluster_density`**: Density weighted by cluster size

## Interpretation

- **Density = 1.0**: On average, each gene connects to 1 term
- **Density = 2.0**: On average, each gene connects to 2 terms
- **Density = 0.5**: On average, each gene connects to 0.5 terms (some genes have connections, some don't)
- **Density = 0.0**: No connections (shouldn't happen if it's a valid cluster)

## Advantages of This Approach

1. **Simple**: Easy to understand and compute
2. **Intuitive**: Directly measures average gene connectivity
3. **Interpretable**: Clear meaning (average connections per gene)
4. **Comparable**: Can be compared across clusters (though cluster size may affect connectivity patterns)

## Why Not Global Network Density?

The old formula `network_density = n_edges / (n_genes × n_terms)` was **incorrect** because:
- It assumes any gene can connect to any term
- But terms have fixed gene sets from libraries
- Not all gene-term pairs are possible

Cluster-level density is correct because:
- It accounts for the library constraint
- It measures connectivity within biologically meaningful groups (clusters)
- It reflects how well genes connect terms within a functional group

