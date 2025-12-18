# Handling Gene List Sizes in P-Value Correction

## Problem

Permutation data was generated for **specific gene list sizes**: 50, 100, 150, 200, ..., 1000 (increments of 50).

However, **users can have any gene list size** between 50-1000 (e.g., 73 genes, 247 genes, 892 genes).

## Solution: Nearest Neighbor Approach

We use a **nearest neighbor** strategy to map any gene list size to the closest available size in the permutation data.

### How It Works

1. **Identify available sizes** from permutation data: [50, 100, 150, 200, ..., 1000]

2. **For any user gene list size**, find the nearest available size:
   - 73 genes → maps to **50** (nearest)
   - 247 genes → maps to **250** (nearest)
   - 892 genes → maps to **900** (nearest)
   - 75 genes → maps to **100** (nearest, since 75 is closer to 100 than 50)

3. **Use the mapped size** in stratification to find the correct null distribution

### Example

```python
# User has 73 genes
gene_list_size = 73

# Find nearest available size
available_sizes = [50, 100, 150, ...]
nearest = 50  # |73 - 50| = 23, |73 - 100| = 27, so 50 is closer

# Use size 50 for correction
stratum = (library, iteration, term_size, overlap, 50)
```

## Implementation Details

### Updated Stratification

Now includes **5 dimensions**:
1. Library (GO BP, GO CC, GO MF, KEGG, Reactome)
2. Iteration bin (1, 2-3, 4-6, 7-10, 11+)
3. Term size bin (<50, 50-100, 100-200, 200-300, 300+)
4. Overlap size bin (3, 4, 5, 6-7, 8+)
5. **Gene list size** (mapped to nearest: 50, 100, 150, ..., 1000)

### Impact on Strata

- **Before**: ~527 strata (without gene list size)
- **After**: ~2,753 strata (with gene list size)
- **Total possible**: ~6,459 strata (some are sparse)

The increase in strata is expected and beneficial because:
- We're accounting for an additional dependency (gene list size correlation: -0.08)
- More precise correction for different gene list sizes
- Sparse strata use fallback to parent strata (without gene list size)

## Why This Approach?

### Alternative Approaches Considered

1. **Binning into ranges** (e.g., 50-100, 100-200):
   - ❌ Less precise
   - ❌ Might mix different behaviors within bins

2. **Interpolation** between sizes:
   - ❌ Complex to implement
   - ❌ Assumes smooth behavior (may not be true)
   - ❌ Requires parametric assumptions

3. **Nearest neighbor** (chosen):
   - ✅ Simple and robust
   - ✅ Uses actual permutation data (no assumptions)
   - ✅ Works well since sizes are close (max distance: 25 genes)

### Validation

The nearest neighbor approach is reasonable because:
- **Maximum distance**: At most 25 genes away (e.g., 75 → 50 or 100)
- **Small effect size**: Gene list size correlation is only -0.08 (moderate)
- **Smooth behavior**: P-value distributions likely change smoothly with gene list size

## Usage

### In Code

```python
from code.pvalue_corrector import load_corrector_from_permutations

corrector = load_corrector_from_permutations("permutation_results/merged_permutation_results.tsv")

# User has 73 genes (not in permutation data)
p_corrected = corrector.correct_pvalue(
    p_raw=0.01,
    library="GO BP",
    iteration=2,
    term_size=150,
    overlap_size=4,
    gene_list_size=73  # Will automatically map to 50
)
```

### In DataFrame Correction

```python
# DataFrame must have 'gene_list_size' column
df['gene_list_size'] = len(user_gene_list)  # e.g., 73

df_corrected = corrector.correct_dataframe(df)
# Automatically handles mapping to nearest size
```

## Edge Cases

### Gene List Size < 50

If user has fewer than 50 genes:
- Maps to **50** (smallest available size)
- Warning: Correction may be less accurate for very small lists

### Gene List Size > 1000

If user has more than 1000 genes:
- Maps to **1000** (largest available size)
- Warning: Correction may be less accurate for very large lists

### Missing Gene List Size

If `gene_list_size` is not provided:
- Uses parent stratum (without gene list size dimension)
- Still works, but less precise

## Recommendations

1. **Always provide gene_list_size** when correcting p-values
2. **For gene lists < 50 or > 1000**: Consider warning users that correction may be less accurate
3. **For future work**: Consider generating permutation data for more sizes (e.g., every 25 genes) for better precision

## Summary

✅ **Problem solved**: Gene list sizes not in permutation data are handled by nearest neighbor mapping

✅ **Implementation**: Updated `pvalue_corrector.py` to include gene list size in stratification

✅ **Impact**: More precise correction accounting for gene list size dependency (correlation: -0.08)

✅ **Robustness**: Fallback to parent strata if mapped size is sparse
