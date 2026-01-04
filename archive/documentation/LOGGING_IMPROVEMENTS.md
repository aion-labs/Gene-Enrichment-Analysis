# Logging Improvements for Permutation Script

## Problem
The permutation script was generating extremely verbose terminal output, making it difficult to monitor progress. Every permutation was printing:
- Library loading messages
- Gene set creation messages  
- Enrichment calculation details
- Iteration details
- Background filtering messages

This resulted in thousands of log lines per second, making it impossible to see actual progress.

## Solution
Modified `code/generate_permutation_distribution.py` to implement quieter logging:

### Changes Made

1. **Separate File and Console Logging**
   - **File handler**: Still logs everything at INFO level to `permutation_distribution.log` (for debugging)
   - **Console handler**: Only shows WARNING and above by default

2. **Progress Messages Still Visible**
   - Added a separate console handler specifically for the main script's logger
   - This handler shows INFO messages from the permutation script (progress updates)
   - Uses a filter to only show messages from `generate_permutation_distribution` module

3. **Suppressed Child Module Verbosity**
   - Child modules (`gene_set_library`, `gene_set`, `enrichment`, `iter_enrichment`, `background_gene_set`) still log to file
   - But their INFO messages are suppressed on console (only WARNING+ shown)

4. **Worker Process Logging**
   - Worker processes (in `run_single_permutation`) also suppress verbose logging
   - Ensures quiet operation even with 32 parallel jobs

## Result

**Before**: Thousands of log lines per second
```
2025-12-09 21:12:11,093 - INFO - Gene Set Library
Reactome (Homo Sapiens)
        1787 terms
        11369 genes
2025-12-09 21:12:11,096 - INFO - Gene Set Library
Reactome (Homo Sapiens)
...
```

**After**: Clean progress messages only
```
Progress: 50/1000 (5.0%) | Elapsed: 0:05:23 | Rate: 0.16 perm/s | ETA: 1:38:45 | Success: 50, Failed: 0
Progress: 100/1000 (10.0%) | Elapsed: 0:10:45 | Rate: 0.16 perm/s | ETA: 1:33:20 | Success: 100, Failed: 0
```

## What You'll See

### Console Output (Terminal)
- ✅ Progress messages (every 1% or 50 permutations)
- ✅ Warnings and errors
- ❌ Verbose library/gene set/enrichment details

### Log File (`permutation_distribution.log`)
- ✅ Everything (full details for debugging)
- ✅ All INFO messages from all modules
- ✅ Complete trace of what happened

## Benefits

1. **Easy Monitoring**: Can now see actual progress without scrolling through thousands of lines
2. **Still Debuggable**: Full logs still written to file if issues occur
3. **No Information Loss**: All details preserved in log file
4. **Better Performance**: Less I/O overhead from excessive console writes

## Usage

No changes needed - the script automatically uses quieter logging. The next time you run:

```bash
python code/generate_permutation_distribution.py --n-permutations 1000 --n-jobs 32
```

You'll see clean progress output instead of verbose logs.








