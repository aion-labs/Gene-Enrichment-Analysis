# Regenerating Permutation Data

## Issue Identified

The permutation data was generated in **multiple separate runs** and then combined. This is **invalid** because:

1. **Each run generates its own random gene lists** - Even with seeds, separate runs may produce different random gene lists for the same permutation number
2. **Different runs may have used different:**
   - Background gene sets
   - Library versions
   - Random number generator states
   - Python versions

3. **Combining permutations from different runs mixes data from different underlying gene sets**, invalidating the statistical analysis

## Solution: Regenerate All Permutations in a Single Run

All permutations must be generated in **one consistent run** to ensure:
- All permutations use the same background gene set
- All permutations use the same library versions
- All permutations use the same random number generator state
- Permutation N always generates the same random gene list (due to seed=perm_idx)

## What Needs to Be Regenerated

### Current Data (INVALID - DO NOT USE):
- `results/permutations_results-50-500/` (sizes 50-500)
- `results/permutation_results-550-1000/` (sizes 550-1000)
- `results/permutation_results-FirstRun-50-to-1000/` (all sizes 50-1000)
- `results/combined_permutation_data/` (combined data - INVALID)
- `results/permutation_cluster_statistics.tsv` (cluster stats - INVALID)
- `results/permutation_cluster_statistics_parquet/` (parquet files - INVALID)

### Libraries to Include (13 total):

**Complete coverage (100%):**
1. C2: CGP: Chemical & genetic perturbations
2. C2: CP: Canonical pathways
3. C2: CP: Pathway Interaction Database
4. C2: CP: WikiPathways
5. GO BP
6. GO CC
7. GO MF
8. H: Hallmark Gene Sets
9. KEGG
10. Protein Interaction
11. Reactome

**High partial coverage (85-90%):**
12. C2: CP: BioCarta (85% - missing sizes 50, 100, 150)
13. C2: CP: KEGG MEDICUS (90% - missing sizes 50, 150)

### Gene List Sizes:
All sizes from 50 to 1000 in steps of 50: 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000

### Permutations per Size:
1000 permutations per size (20,000 total permutations)

## Regeneration Steps

### Step 1: Fix Column Filtering Issue

**IMPORTANT:** The script has been updated to keep **all columns** (not just essential ones) to match the format of the 50-500 run. The previous 550-1000 run was missing columns like:
- `iteration -log(p-value)`
- `Description`
- `iteration overlapping genes`
- `Full list overlapping genes`
- `Full list p-value`
- `Regular FDR`
- `Full list overlapping genes (gene names)`

The script now keeps all columns for consistency. This has already been fixed in `code/generate_permutation_distribution.py`.

### Step 2: Update Library Configuration

The script `code/generate_permutation_distribution.py` already includes most libraries, but **library names must match exactly** what's expected. Check that the library names in the `LIBRARIES` dictionary match:

**Required library names (must match exactly):**
- "H: Hallmark Gene Sets" ✓
- "C2: CGP: Chemical & genetic perturbations" ✓
- "C2: CP: Canonical pathways" ✓
- "C2: CP: BioCarta" ✓
- "C2: CP: KEGG MEDICUS" ✓
- "C2: CP: Pathway Interaction Database" ✓
- "C2: CP: WikiPathways" ✓
- "Protein Interaction" ✓

**Library name mapping needed:**

The script currently has different names than what was used in previous runs. Update the `LIBRARIES` dictionary in `code/generate_permutation_distribution.py`:

**Current (in script) → Required (from previous runs):**
- "C2: CP: Reactome Pathways" → "Reactome"
- "C2: CP: KEGG Legacy" → "KEGG"
- "C5: GO: Biological Process" → "GO BP"
- "C5: GO: Cellular Component" → "GO CC"
- "C5: GO: Molecular Function" → "GO MF" (already correct in some runs)

**Action:** Update the `LIBRARIES` dictionary to use these exact names:
```python
LIBRARIES = {
    "H: Hallmark Gene Sets": "h.all.v2025.1.Hs.symbols.gmt",
    "C2: CGP: Chemical & genetic perturbations": "c2.cgp.v2025.1.Hs.symbols.gmt",
    "C2: CP: Canonical pathways": "c2.cp.v2025.1.Hs.symbols.gmt",
    "C2: CP: BioCarta": "c2.cp.biocarta.v2025.1.Hs.symbols.gmt",
    "C2: CP: KEGG MEDICUS": "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt",
    "C2: CP: Pathway Interaction Database": "c2.cp.pid.v2025.1.Hs.symbols.gmt",
    "C2: CP: WikiPathways": "c2.cp.wikipathways.v2025.1.Hs.symbols.gmt",
    "Reactome": "c2.cp.reactome.v2025.1.Hs.symbols.gmt",  # Changed from "C2: CP: Reactome Pathways"
    "KEGG": "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",  # Changed from "C2: CP: KEGG Legacy"
    "GO BP": "c5.go.bp.v2025.1.Hs.symbols.gmt",  # Changed from "C5: GO: Biological Process"
    "GO CC": "c5.go.cc.v2025.1.Hs.symbols.gmt",  # Changed from "C5: GO: Cellular Component"
    "GO MF": "c5.go.mf.v2025.1.Hs.symbols.gmt",  # Changed from "C5: GO: Molecular Function"
    "Protein Interaction": "stringdb_protein_interactions.gmt",
}
```

### Step 3: Run Permutation Generation (Single Run)

```bash
# Activate virtual environment
source venv/bin/activate  # or .venv/bin/activate

# Run permutation generation for all sizes in one run
python code/generate_permutation_distribution.py \
    --n-permutations 1000 \
    --n-jobs 32 \
    --resume

# This will:
# - Generate all 20,000 permutations (1000 per size × 20 sizes)
# - Use consistent seeds (perm_idx) for reproducibility
# - Save to results/permutation_results/ (or configured output directory)
```

### Step 4: Verify Consistency

After generation, verify that:
- All permutations for a given size use the same background
- Permutation N always generates the same random gene list (test by regenerating a few)
- All libraries are included in the results

### Step 5: Extract Cluster Statistics

```bash
# Extract cluster statistics from the single run
python extract_cluster_stats_from_combined_data.py \
    --combined-data-dir results/permutation_results \
    --output results/permutation_cluster_statistics.tsv
```

### Step 6: Create Parquet Files

```bash
# Create parquet files for benchmarking
python prepare_parquet_cluster_stats.py \
    --input results/permutation_cluster_statistics.tsv \
    --output-dir results/permutation_cluster_statistics_parquet
```

## Important Notes

1. **Do NOT combine data from multiple runs** - Always generate all permutations in a single run
2. **Use `--resume` flag** - This allows you to resume if the process is interrupted, but only within the same run
3. **Verify seeds are working** - Test that permutation 1 always generates the same gene list
4. **Check library versions** - Ensure all libraries are the same version throughout the run
5. **Monitor for errors** - Some permutations may fail (no significant results), which is normal

## Estimated Time

- **Permutation generation**: ~20-40 hours depending on hardware (20,000 permutations × ~5-10 seconds each)
- **Cluster statistics extraction**: ~1-2 hours
- **Parquet file creation**: ~5-10 minutes

## Backup Current Data

Before regenerating, consider backing up the current (invalid) data:

```bash
# Backup current data
mkdir -p archive/invalid_permutation_data_$(date +%Y%m%d)
mv results/permutations_results-50-500 archive/invalid_permutation_data_$(date +%Y%m%d)/
mv results/permutation_results-550-1000 archive/invalid_permutation_data_$(date +%Y%m%d)/
mv results/permutation_results-FirstRun-50-to-1000 archive/invalid_permutation_data_$(date +%Y%m%d)/
mv results/combined_permutation_data archive/invalid_permutation_data_$(date +%Y%m%d)/
mv results/permutation_cluster_statistics.tsv archive/invalid_permutation_data_$(date +%Y%m%d)/
mv results/permutation_cluster_statistics_parquet archive/invalid_permutation_data_$(date +%Y%m%d)/
```

