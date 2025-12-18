# P-Value Corrector Integration Guide

## How It Works

The `pvalue_corrector.py` works in **two phases**:

1. **Model Building (One-Time Setup)**: Load permutation data and build CDFs for each stratum
2. **Correction (During iGEA)**: Use the model to correct p-values as they're computed

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│ Phase 1: Model Building (One-Time)                          │
│                                                              │
│  Permutation Data → Load → Build CDFs → Save Model          │
│  (1.25M rows)      (TSV)   (2,753 strata)  (Optional)       │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ Phase 2: Correction (During iGEA)                           │
│                                                              │
│  iGEA computes p-value → Corrector corrects it → Use result│
│  (raw p-value)      (model)    (corrected p-value)          │
└─────────────────────────────────────────────────────────────┘
```

## Two Integration Approaches

### Approach 1: Post-Processing (Easier, Recommended Initially)

**How it works:**
- Run iGEA normally (computes raw p-values)
- After iGEA completes, apply correction to all results
- Add corrected p-values as a new column

**Advantages:**
- ✅ No changes to core iGEA logic
- ✅ Easy to implement
- ✅ Can be applied retroactively to existing results
- ✅ Can compare raw vs. corrected

**Disadvantages:**
- ❌ Two-step process
- ❌ Can't use corrected p-values for iteration decisions

### Approach 2: Integrated (More Complex, Better Long-Term)

**How it works:**
- Initialize corrector once at startup
- Pass corrector to `IterativeEnrichment`
- Correct p-values immediately after computation
- Use corrected p-values for iteration decisions

**Advantages:**
- ✅ Single-step process
- ✅ Can use corrected p-values for stopping criteria
- ✅ More accurate iteration decisions

**Disadvantages:**
- ❌ Requires modifying `IterativeEnrichment` class
- ❌ More complex implementation

---

## Implementation: Approach 1 (Post-Processing)

### Step 1: Initialize Corrector (One-Time)

```python
from code.pvalue_corrector import load_corrector_from_permutations
from pathlib import Path

# Initialize corrector (do this once, can reuse)
project_root = Path(__file__).parent.parent
permutation_file = project_root / "permutation_results" / "merged_permutation_results.tsv"

corrector = load_corrector_from_permutations(str(permutation_file))
# This loads 1.25M rows and builds 2,753 CDFs (takes ~10-30 seconds)
```

### Step 2: Run iGEA Normally

```python
from code.iter_enrichment import IterativeEnrichment

# Run iGEA as usual
iter_enrich = IterativeEnrichment(
    gene_set=gene_set,
    gene_set_library=library,
    background_gene_set=background,
    # ... other parameters
)

# Get results (with raw p-values)
results = iter_enrich.results
```

### Step 3: Apply Correction

```python
import pandas as pd

# Convert to DataFrame
df = iter_enrich.to_dataframe()

# Add gene list size (needed for correction)
df['gene_list_size'] = len(gene_set.genes)

# Apply correction
df_corrected = corrector.correct_dataframe(df)

# Now df_corrected has 'corrected p-value' column
```

### Complete Example

```python
#!/usr/bin/env python3
"""Example: Post-processing correction"""

from pathlib import Path
from code.pvalue_corrector import load_corrector_from_permutations
from code.iter_enrichment import IterativeEnrichment
from code.gene_set import GeneSet
from code.gene_set_library import GeneSetLibrary
from code.background_gene_set import BackgroundGeneSet

# 1. Initialize corrector (one-time, can be reused)
project_root = Path(__file__).parent.parent
permutation_file = project_root / "permutation_results" / "merged_permutation_results.tsv"
corrector = load_corrector_from_permutations(str(permutation_file))

# 2. Setup iGEA
gene_set = GeneSet(["GENE1", "GENE2", ...], background_genes)
library = GeneSetLibrary.load("path/to/library.gmt")
background = BackgroundGeneSet.load("path/to/background.txt")

# 3. Run iGEA
iter_enrich = IterativeEnrichment(
    gene_set=gene_set,
    gene_set_library=library,
    background_gene_set=background,
    p_threshold=0.01,
    max_iterations=10,
    # ... other params
)

# 4. Get results and apply correction
df = iter_enrich.to_dataframe()
df['gene_list_size'] = len(gene_set.genes)
df_corrected = corrector.correct_dataframe(df)

# 5. Save corrected results
df_corrected.to_csv("results_with_correction.tsv", sep="\t", index=False)

# 6. Compare raw vs. corrected
print("Raw p-values:", df['iteration p-value'].head())
print("Corrected p-values:", df_corrected['corrected p-value'].head())
```

---

## Implementation: Approach 2 (Integrated)

### Step 1: Modify `IterativeEnrichment.__init__`

Add optional `pvalue_corrector` parameter:

```python
class IterativeEnrichment:
    def __init__(
        self,
        gene_set: GeneSet,
        gene_set_library: GeneSetLibrary,
        background_gene_set: BackgroundGeneSet,
        # ... existing parameters ...
        pvalue_corrector: Optional[iGEAPValueCorrector] = None,  # NEW
    ):
        # ... existing initialization ...
        self.pvalue_corrector = pvalue_corrector
        self.gene_list_size = len(gene_set.genes)  # Store for correction
```

### Step 2: Modify `_compute_enrichment` to Correct P-Values

```python
def _compute_enrichment(self) -> List[Dict[str, Any]]:
    # ... existing code ...
    
    # After getting top result
    top = filtered_results[0]
    pval = top.get("p-value")
    
    # CORRECT P-VALUE HERE
    if self.pvalue_corrector is not None:
        # Extract needed information
        overlap_size_str = top.get("overlap_size", "0/0")
        overlap_count = int(overlap_size_str.split("/")[0])
        term_size = int(overlap_size_str.split("/")[1])
        
        # Correct p-value
        pval_corrected = self.pvalue_corrector.correct_pvalue(
            p_raw=pval,
            library=self.gene_set_library.name,
            iteration=iteration,
            term_size=term_size,
            overlap_size=overlap_count,
            gene_list_size=self.gene_list_size
        )
        
        # Use corrected p-value for threshold check
        if pval_corrected >= self.p_threshold:
            break
        
        # Store both raw and corrected
        pval = pval_corrected
    else:
        # Original behavior (no correction)
        if pval >= self.p_threshold:
            break
    
    # ... rest of code ...
    
    record = {
        "Iteration": iteration,
        "Term": format_term_name(top.get("term", "")),
        # ... other fields ...
        "iteration p-value": pval,  # Now corrected if corrector provided
        "raw p-value": top.get("p-value"),  # Keep original
        # ... rest ...
    }
```

### Step 3: Use Corrected P-Values in Results

```python
# In to_dataframe() or similar methods
if self.pvalue_corrector is not None:
    df['raw p-value'] = df['iteration p-value']  # Keep original
    df['iteration p-value'] = df['corrected p-value']  # Use corrected
```

### Complete Example

```python
#!/usr/bin/env python3
"""Example: Integrated correction"""

from code.pvalue_corrector import load_corrector_from_permutations
from code.iter_enrichment import IterativeEnrichment

# 1. Initialize corrector
corrector = load_corrector_from_permutations("permutation_results/merged_permutation_results.tsv")

# 2. Run iGEA with corrector
iter_enrich = IterativeEnrichment(
    gene_set=gene_set,
    gene_set_library=library,
    background_gene_set=background,
    pvalue_corrector=corrector,  # Pass corrector
    # ... other params
)

# Results now have corrected p-values
df = iter_enrich.to_dataframe()
# df['iteration p-value'] is already corrected
```

---

## Where to Initialize Corrector

### Option A: Global Singleton (Recommended)

Create a module-level corrector that's initialized once:

```python
# code/pvalue_corrector_singleton.py
from code.pvalue_corrector import load_corrector_from_permutations
from pathlib import Path

_corrector = None

def get_corrector():
    """Get or initialize the global corrector."""
    global _corrector
    if _corrector is None:
        project_root = Path(__file__).parent.parent
        permutation_file = project_root / "permutation_results" / "merged_permutation_results.tsv"
        _corrector = load_corrector_from_permutations(str(permutation_file))
    return _corrector
```

**Usage:**
```python
from code.pvalue_corrector_singleton import get_corrector

corrector = get_corrector()  # Initialized once, reused
```

### Option B: Initialize in Streamlit App

```python
# In streamlit_app.py
import streamlit as st
from code.pvalue_corrector import load_corrector_from_permutations

@st.cache_resource  # Cache so it's only loaded once
def load_pvalue_corrector():
    permutation_file = ROOT / "permutation_results" / "merged_permutation_results.tsv"
    return load_corrector_from_permutations(str(permutation_file))

# In main()
corrector = load_pvalue_corrector()
```

### Option C: Lazy Initialization

Initialize only when needed:

```python
class IterativeEnrichment:
    def __init__(self, ..., use_pvalue_correction: bool = False):
        self.use_pvalue_correction = use_pvalue_correction
        self._corrector = None
    
    @property
    def corrector(self):
        if self.use_pvalue_correction and self._corrector is None:
            self._corrector = load_corrector_from_permutations(...)
        return self._corrector
```

---

## Performance Considerations

### Model Building (One-Time)
- **Time**: ~10-30 seconds (depends on hardware)
- **Memory**: ~500MB-1GB (stores permutation data + CDFs)
- **When**: Once at startup or first use

### Correction (Per P-Value)
- **Time**: ~0.001-0.01ms per p-value (very fast)
- **Memory**: Negligible
- **When**: Every time a p-value is computed

### Optimization Tips

1. **Cache the corrector**: Use singleton or Streamlit `@st.cache_resource`
2. **Lazy loading**: Only load when correction is requested
3. **Batch correction**: Use `correct_dataframe()` for multiple p-values at once

---

## Recommended Workflow

### For Development/Testing

1. **Start with Approach 1** (post-processing)
   - Easier to implement
   - Can test without modifying core code
   - Can compare raw vs. corrected

2. **Test on sample data**
   - Verify corrections are reasonable
   - Check that corrected p-values are more uniform

3. **Move to Approach 2** (integrated) if needed
   - Only if you want corrected p-values to affect iteration decisions
   - Requires more testing

### For Production

1. **Initialize corrector once** at application startup
2. **Use Approach 1** initially (safer, easier to debug)
3. **Consider Approach 2** for better accuracy (if needed)

---

## Example: Adding to Streamlit App

```python
# In streamlit_app.py

import streamlit as st
from code.pvalue_corrector import load_corrector_from_permutations
from pathlib import Path

ROOT = Path(__file__).parent.parent

@st.cache_resource
def get_pvalue_corrector():
    """Load p-value corrector (cached)."""
    permutation_file = ROOT / "permutation_results" / "merged_permutation_results.tsv"
    if permutation_file.exists():
        return load_corrector_from_permutations(str(permutation_file))
    return None

# In main(), after running iGEA
if st.checkbox("Apply p-value correction"):
    corrector = get_pvalue_corrector()
    if corrector:
        # Apply correction to all iterative results
        for lib_name, iter_enrich in state.iter_enrich.items():
            df = iter_enrich.to_dataframe()
            df['gene_list_size'] = len(state.gene_set.genes)
            df_corrected = corrector.correct_dataframe(df)
            
            # Update results with corrected p-values
            # (implementation depends on how you want to display/store)
            st.dataframe(df_corrected[['Term', 'iteration p-value', 'corrected p-value']])
    else:
        st.warning("P-value corrector not available (permutation file missing)")
```

---

## Summary

- **Model Building**: One-time setup (load permutation data, build CDFs)
- **Correction**: Apply during or after iGEA
- **Approach 1 (Post-Processing)**: Easier, recommended initially
- **Approach 2 (Integrated)**: More complex, better long-term
- **Performance**: Fast correction (~0.001ms per p-value), model building takes 10-30s

Choose the approach that fits your needs!
