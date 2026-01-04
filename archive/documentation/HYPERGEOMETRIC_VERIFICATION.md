# Hypergeometric P-Value Calculation Verification

## Test Setup
- **Input Gene Set**: Hypoxia genes (67 valid genes from 75 total)
- **Library**: Reactome Pathways (1,787 terms, 11,369 unique genes)
- **Background**: All genes (23,391 genes)
- **Library-Specific Background**: 11,350 genes (intersection of background and Reactome genes)

## Hypergeometric Parameters

The hypergeometric test uses these parameters:
- **M**: Total population size (library_background_size = 11,350)
- **n**: Number of successes in population (term size)
- **N**: Sample size (gene set size = 67)
- **k**: Number of successes in sample (overlap size)

## Example Calculations

### 1. REACTOME_METABOLISM_OF_CARBOHYDRATES_AND_CARBOHYDRATE_DERIVATIVES
- **Term size (n)**: 288
- **Overlap (k)**: 5 genes (ALDOC, GAPDH, GYS1, PGAM1, PGK1)
- **Parameters**: M=11,350, n=288, N=67, k=5
- **Scipy call**: `hypergeom.sf(4, 11350, 288, 67)`
- **P-value**: 2.73e-02
- **Interpretation**: P(X ≥ 5) where X ~ Hypergeometric(M=11,350, n=288, N=67)

### 2. REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53
- **Term size (n)**: 360
- **Overlap (k)**: 5 genes (CCNB1, CYCS, DDIT4, NDRG1, RRAGD)
- **Parameters**: M=11,350, n=360, N=67, k=5
- **Scipy call**: `hypergeom.sf(4, 11350, 360, 67)`
- **P-value**: 6.10e-02
- **Interpretation**: P(X ≥ 5) where X ~ Hypergeometric(M=11,350, n=360, N=67)

### 3. REACTOME_COLLAGEN_FORMATION
- **Term size (n)**: 90
- **Overlap (k)**: 4 genes (LOX, P4HA1, P4HA2, PLOD1)
- **Parameters**: M=11,350, n=90, N=67, k=4
- **Scipy call**: `hypergeom.sf(3, 11350, 90, 67)`
- **P-value**: 1.94e-03
- **Interpretation**: P(X ≥ 4) where X ~ Hypergeometric(M=11,350, n=90, N=67)

## Key Verification Points

### ✅ Library-Specific Background Size
- **Correct**: Uses intersection of background genes and library genes (11,350)
- **Not**: Full background size (23,391) or library size (11,369)

### ✅ Survival Function Usage
- **Correct**: `hypergeom.sf(k-1, M, n, N)` gives P(X ≥ k)
- **Not**: `hypergeom.sf(k, M, n, N)` which would give P(X > k)

### ✅ Parameter Order
- **Correct**: `hypergeom.sf(k, M, n, N)` where:
  - k = overlap - 1
  - M = library_background_size
  - n = term_size
  - N = gene_set_size

### ✅ Biological Interpretation
The p-values correctly reflect the enrichment of hypoxia genes in metabolic and regulatory pathways, with the most significant enrichment in collagen formation (p = 1.94e-03).

## Formula Verification

The hypergeometric distribution models drawing k successes from a population of size M containing n successes, when drawing N items:

**P(X = k) = C(n,k) × C(M-n, N-k) / C(M,N)**

For enrichment analysis, we want **P(X ≥ k)**, calculated using the survival function:
**hypergeom.sf(k-1, M, n, N)**

This gives the probability of observing k or more successes by chance, which is exactly what we want for gene set enrichment analysis.

## Conclusion

The hypergeometric p-value calculation is implemented correctly in the codebase. The key improvements include:

1. **Library-specific background size**: Uses intersection of background and library genes
2. **Correct survival function**: Uses `k-1` to get P(X ≥ k)
3. **Proper parameter order**: Follows scipy.stats.hypergeom.sf convention
4. **Biological relevance**: Results make sense for hypoxia gene enrichment

The implementation in `code/enrichment.py` line 78-80 is mathematically correct:

```python
p_value = hypergeom.sf(
    n_overlap - 1, library_background_size, n_term_genes, gene_set.size
)
```
