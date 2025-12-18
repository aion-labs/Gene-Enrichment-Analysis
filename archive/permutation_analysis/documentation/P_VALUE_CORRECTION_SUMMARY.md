# P-Value Correction for iGEA - Summary

## Problem

iGEA p-values are dependent on:
1. **Library** (strongest effect)
2. **Iteration number** (correlation: 0.33)
3. **Term size** (correlation: -0.11)
4. **Overlap size** (weak correlation: -0.035)

This violates the assumption of uniform p-values under the null hypothesis, making standard multiple testing corrections (like FDR) inappropriate.

## Solution: Stratified Empirical Null Distribution

### Approach

Create separate null distributions for each combination of:
- Library (5 levels)
- Iteration bin (1, 2-3, 4-6, 7-10, 11+)
- Term size bin (<50, 50-100, 100-200, 200-300, 300+)
- Overlap size bin (3, 4, 5, 6-7, 8+)

For each stratum, build an empirical CDF from permutation data, then transform observed p-values to be uniform under the null.

### Implementation

See `code/pvalue_corrector.py` for the full implementation.

**Key features:**
- Non-parametric (no distributional assumptions)
- Handles sparse strata with fallback to parent strata
- Transforms p-values to uniform scale under null
- Easy to use with existing iGEA results

### Usage

```python
from code.pvalue_corrector import load_corrector_from_permutations

# Load corrector from permutation results
corrector = load_corrector_from_permutations(
    "permutation_results/merged_permutation_results.tsv"
)

# Correct a single p-value
p_corrected = corrector.correct_pvalue(
    p_raw=0.01,
    library="GO BP",
    iteration=2,
    term_size=150,
    overlap_size=4
)

# Or correct an entire DataFrame
df_corrected = corrector.correct_dataframe(df)
```

### Results

From 1.25M permutation results:
- Built **527 stratum CDFs** (out of 600 possible strata)
- 73 strata were too sparse (<100 samples) and use fallback correction
- Correction accounts for all major dependencies

### Interpretation

The corrected p-value represents the percentile rank of the raw p-value within its stratum's null distribution. A corrected p-value of 0.05 means the raw p-value is at the 5th percentile of the null distribution for that specific combination of library, iteration, term size, and overlap.

### Next Steps

1. **Integrate into iGEA**: Add corrected p-values to output
2. **Validate**: Check that corrected p-values are uniform under null
3. **Compare**: Test on real data to ensure power is maintained
4. **Refine**: Adjust binning strategy if needed

### Files Created

- `P_VALUE_CORRECTION_RECOMMENDATIONS.md`: Detailed recommendations and alternative approaches
- `code/pvalue_corrector.py`: Implementation of stratified correction
- `P_VALUE_CORRECTION_SUMMARY.md`: This summary
