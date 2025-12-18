# P-Value Correction Recommendations for iGEA

## Executive Summary

Based on analysis of permutation results, iGEA p-values show significant dependencies on:
1. **Library** (strongest effect, different distributions per library)
2. **Iteration number** (correlation: 0.33, p-values increase with iteration)
3. **Term size** (correlation: -0.11, larger terms have lower p-values)
4. **Overlap size** (weak correlation: -0.035)
5. **Gene list size** (correlation: -0.08, included in stratification)

**Note**: Permutation data has gene list sizes at 50, 100, 150, ..., 1000. For intermediate sizes, we use nearest neighbor mapping (see `GENE_LIST_SIZE_HANDLING.md`).

## Recommended Approaches

### Approach 1: Stratified Empirical Null Distribution (RECOMMENDED)

**Concept**: Create separate null distributions for each combination of covariates, then use empirical quantiles to correct p-values.

**Implementation**:
1. **Stratify permutation data** by:
   - Library (5 levels: GO BP, GO CC, GO MF, KEGG, Reactome)
   - Iteration number (binned: 1, 2-3, 4-6, 7-10, 11+)
   - Term size (binned: <50, 50-100, 100-200, 200-300, 300+)
   - Overlap size (binned: 3, 4, 5, 6-7, 8+)
   - Gene list size (mapped to nearest: 50, 100, 150, ..., 1000)

2. **For each stratum**, compute:
   - Empirical CDF of p-values
   - Quantiles (e.g., 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)

3. **Correct observed p-values**:
   - Find the appropriate stratum for each test
   - Use empirical quantile function to convert raw p-value to corrected p-value
   - Formula: `p_corrected = empirical_quantile(p_raw, stratum)`

**Advantages**:
- Non-parametric, makes no distributional assumptions
- Accounts for all dependencies
- Directly uses permutation data
- Easy to interpret

**Disadvantages**:
- Requires sufficient data in each stratum
- May need smoothing/interpolation for sparse strata

**Code Structure**:
```python
class StratifiedNullDistribution:
    def __init__(self, permutation_df):
        # Create strata
        # Build empirical CDFs per stratum
        
    def correct_pvalue(self, p_raw, library, iteration, term_size, overlap_size):
        # Find stratum
        # Get empirical quantile
        # Return corrected p-value
```

---

### Approach 2: Parametric Beta Distribution Modeling

**Concept**: Fit beta distributions to p-values within each stratum, then use the fitted distribution for correction.

**Implementation**:
1. **Stratify** permutation data (same as Approach 1)

2. **For each stratum**, fit a beta distribution:
   - Beta distribution is natural for p-values (bounded [0,1])
   - Use method of moments or MLE
   - Parameters: α (shape1), β (shape2)

3. **Correct p-values**:
   - For observed p-value, find stratum
   - Use fitted beta CDF to get corrected p-value
   - Optionally transform to uniform scale: `p_corrected = Beta.cdf(p_raw, α, β)`

**Advantages**:
- Smooth, continuous correction
- Can extrapolate beyond observed values
- Statistically principled

**Disadvantages**:
- Assumes beta distribution (may not fit perfectly)
- More complex implementation

**Code Structure**:
```python
from scipy.stats import beta

class BetaNullModel:
    def __init__(self, permutation_df):
        # Fit beta distributions per stratum
        self.models = {}  # {stratum: (alpha, beta)}
        
    def correct_pvalue(self, p_raw, library, iteration, term_size, overlap_size):
        stratum = self._get_stratum(library, iteration, term_size, overlap_size)
        alpha, beta = self.models[stratum]
        return beta.cdf(p_raw, alpha, beta)
```

---

### Approach 3: Regression-Based Correction

**Concept**: Model p-values as a function of covariates using regression, then use residuals for correction.

**Implementation**:
1. **Fit regression model** to permutation p-values:
   ```
   logit(p) ~ Library + Iteration + log(TermSize) + OverlapSize + GeneListSize
   ```
   - Use logit link (or log link) since p-values are bounded
   - Include interaction terms if significant

2. **Extract residuals**:
   - Residuals represent deviation from expected null p-value
   - Transform residuals to corrected p-values

3. **Alternative**: Use quantile regression to model different quantiles

**Advantages**:
- Captures interactions between covariates
- Smooth, continuous correction
- Can handle continuous covariates directly

**Disadvantages**:
- More complex model selection
- Requires assumptions about error distribution

**Code Structure**:
```python
from sklearn.linear_model import QuantileRegressor
import statsmodels.api as sm

class RegressionNullModel:
    def __init__(self, permutation_df):
        # Fit quantile regression models
        # Or logistic regression with logit link
        
    def correct_pvalue(self, p_raw, library, iteration, term_size, overlap_size):
        # Predict expected null p-value
        # Compute residual
        # Transform to corrected p-value
```

---

### Approach 4: Local FDR Estimation

**Concept**: Estimate local false discovery rate (FDR) using the null distribution, then convert to corrected p-values.

**Implementation**:
1. **Estimate null distribution** (using Approach 1 or 2)

2. **For each observed p-value**:
   - Compute density under null: `f_null(p)`
   - Estimate mixture model: `f_observed(p) = π₀ * f_null(p) + (1-π₀) * f_alt(p)`
   - Compute local FDR: `lfdr(p) = π₀ * f_null(p) / f_observed(p)`

3. **Convert lFDR to corrected p-value**:
   - Use relationship: `p_corrected ≈ lFDR(p)`

**Advantages**:
- Directly estimates false discovery rate
- Can handle mixture of null and alternative
- Well-established statistical framework

**Disadvantages**:
- Requires estimating mixture proportions
- More computationally intensive

---

### Approach 5: Hybrid: Stratified + Smoothing

**Concept**: Combine stratified approach with smoothing to handle sparse strata.

**Implementation**:
1. **Create strata** (as in Approach 1)

2. **For sparse strata**, use:
   - Kernel smoothing across similar strata
   - Hierarchical modeling (shrinkage toward parent strata)
   - Example: If stratum (Library=GO BP, Iter=5, Term=100-200, Overlap=4) is sparse,
     borrow information from:
     - Same library, nearby iterations
     - Same library, similar term sizes
     - Similar libraries

3. **Apply correction** using smoothed empirical CDF

**Advantages**:
- Handles data sparsity
- More robust than pure stratification
- Still non-parametric

**Disadvantages**:
- More complex implementation
- Requires tuning smoothing parameters

---

## Recommended Implementation Strategy

### Phase 1: Start with Approach 1 (Stratified Empirical)

1. **Create stratification function**:
   ```python
   def create_stratum(library, iteration, term_size, overlap_size):
       iter_bin = bin_iteration(iteration)
       term_bin = bin_term_size(term_size)
       overlap_bin = bin_overlap_size(overlap_size)
       return (library, iter_bin, term_bin, overlap_bin)
   ```

2. **Build empirical CDFs**:
   ```python
   from scipy.stats import rv_histogram
   
   for stratum in all_strata:
       stratum_data = permutation_df[permutation_df['stratum'] == stratum]['p_value']
       # Create histogram
       # Convert to CDF
       stratum_cdf = rv_histogram(histogram)
   ```

3. **Correction function**:
   ```python
   def correct_pvalue(p_raw, library, iteration, term_size, overlap_size):
       stratum = create_stratum(library, iteration, term_size, overlap_size)
       if stratum in cdf_dict:
           # Use empirical CDF to get corrected p-value
           # Option 1: Percentile rank
           p_corrected = percentile_rank(p_raw, stratum_cdf)
           # Option 2: Uniform transformation
           p_corrected = stratum_cdf.cdf(p_raw)
       else:
           # Fallback: use nearest stratum or global correction
           p_corrected = fallback_correction(p_raw)
       return p_corrected
   ```

### Phase 2: Validate and Refine

1. **Check coverage**: Verify corrected p-values are approximately uniform under null
2. **Check power**: Ensure correction doesn't eliminate true signals
3. **Handle edge cases**: Sparse strata, extreme values

### Phase 3: Optional Enhancement

- Add smoothing (Approach 5) if needed
- Consider parametric models (Approach 2) for smoother correction
- Implement local FDR (Approach 4) for additional insights

---

## Key Considerations

### 1. Binning Strategy

**Iteration**: 
- Fine-grained for early iterations (1, 2, 3, 4, 5)
- Coarser for later iterations (6-10, 11-15, 16+)

**Term Size**:
- Use quantile-based bins to ensure sufficient data per bin
- Consider: <25th percentile, 25-50th, 50-75th, 75-90th, 90th+

**Overlap Size**:
- Keep fine-grained (3, 4, 5, 6, 7, 8+) since distribution is discrete

### 2. Handling Sparse Strata

- Minimum data requirement: At least 100 observations per stratum
- For sparse strata:
  - Use parent stratum (remove one dimension)
  - Use kernel smoothing
  - Use global correction as fallback

### 3. Edge Cases

- **Very small p-values** (< 1e-6): May need special handling
- **Very large p-values** (≈ 0.05): Already at threshold, correction may not be needed
- **Missing covariates**: Use most common/default values

### 4. Validation

- **Q-Q plots**: Compare corrected p-values to uniform distribution
- **Type I error rate**: Should be controlled at nominal level (e.g., 0.05)
- **Power analysis**: Compare corrected vs. uncorrected on known true signals

---

## Example Code Structure

```python
import pandas as pd
import numpy as np
from scipy.stats import rv_histogram
from collections import defaultdict

class iGEAPValueCorrector:
    def __init__(self, permutation_df, min_samples_per_stratum=100):
        self.permutation_df = permutation_df
        self.min_samples = min_samples_per_stratum
        self.stratum_cdfs = {}
        self._build_null_distributions()
    
    def _bin_iteration(self, iteration):
        if iteration == 1:
            return 1
        elif iteration <= 3:
            return 2
        elif iteration <= 6:
            return 3
        elif iteration <= 10:
            return 4
        else:
            return 5
    
    def _bin_term_size(self, term_size):
        if term_size < 50:
            return '<50'
        elif term_size < 100:
            return '50-100'
        elif term_size < 200:
            return '100-200'
        elif term_size < 300:
            return '200-300'
        else:
            return '300+'
    
    def _bin_overlap_size(self, overlap_size):
        if overlap_size <= 3:
            return '3'
        elif overlap_size == 4:
            return '4'
        elif overlap_size == 5:
            return '5'
        elif overlap_size <= 7:
            return '6-7'
        else:
            return '8+'
    
    def _create_stratum(self, library, iteration, term_size, overlap_size):
        return (
            library,
            self._bin_iteration(iteration),
            self._bin_term_size(term_size),
            self._bin_overlap_size(overlap_size)
        )
    
    def _build_null_distributions(self):
        # Parse permutation data
        df = self.permutation_df.copy()
        df['overlap_size'] = df['iteration overlapping genes'].str.split('/').str[0].astype(int)
        df['term_size'] = df['iteration overlapping genes'].str.split('/').str[1].astype(int)
        df['p_value'] = df['iteration p-value'].astype(float)
        
        # Create strata
        df['stratum'] = df.apply(
            lambda row: self._create_stratum(
                row['Library'],
                row['Iteration'],
                row['term_size'],
                row['overlap_size']
            ),
            axis=1
        )
        
        # Build CDF for each stratum
        for stratum, stratum_data in df.groupby('stratum'):
            p_values = stratum_data['p_value'].values
            
            if len(p_values) >= self.min_samples:
                # Create histogram
                counts, bins = np.histogram(p_values, bins=50, range=(0, 0.05))
                # Normalize
                counts = counts / counts.sum()
                # Create CDF
                cdf = rv_histogram((counts, bins))
                self.stratum_cdfs[stratum] = cdf
            else:
                # Mark as sparse, will use fallback
                self.stratum_cdfs[stratum] = None
    
    def correct_pvalue(self, p_raw, library, iteration, term_size, overlap_size):
        """
        Correct a p-value using the stratified null distribution.
        
        Returns:
            corrected_pvalue: P-value corrected for dependencies
        """
        stratum = self._create_stratum(library, iteration, term_size, overlap_size)
        
        if stratum in self.stratum_cdfs and self.stratum_cdfs[stratum] is not None:
            # Use stratum-specific CDF
            cdf = self.stratum_cdfs[stratum]
            # Transform to uniform scale
            p_corrected = cdf.cdf(p_raw)
        else:
            # Fallback: use parent stratum or global correction
            p_corrected = self._fallback_correction(p_raw, stratum)
        
        return p_corrected
    
    def _fallback_correction(self, p_raw, stratum):
        # Try parent strata (remove dimensions one at a time)
        library, iter_bin, term_bin, overlap_bin = stratum
        
        # Try without overlap bin
        parent_stratum = (library, iter_bin, term_bin, None)
        if parent_stratum in self.stratum_cdfs:
            return self.stratum_cdfs[parent_stratum].cdf(p_raw)
        
        # Try without term bin
        parent_stratum = (library, iter_bin, None, None)
        if parent_stratum in self.stratum_cdfs:
            return self.stratum_cdfs[parent_stratum].cdf(p_raw)
        
        # Use library-specific global correction
        library_stratum = (library, None, None, None)
        if library_stratum in self.stratum_cdfs:
            return self.stratum_cdfs[library_stratum].cdf(p_raw)
        
        # Final fallback: no correction
        return p_raw
```

---

## Next Steps

1. **Implement Approach 1** (stratified empirical) as baseline
2. **Validate** on permutation data (should give uniform corrected p-values)
3. **Test** on real data to ensure power is maintained
4. **Refine** binning strategy based on data distribution
5. **Consider** adding smoothing or parametric models if needed

---

## References

- Efron, B. (2010). Large-scale inference: empirical Bayes methods for estimation, testing, and prediction.
- Storey, J. D. (2002). A direct approach to false discovery rates.
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing.
