# Permutation Distribution Script - Usage Guide

## Overview
This script generates a null distribution of p-values for iterative enrichment analysis (iGEA) by running 20,000 permutations (20 gene list sizes × 1000 permutations each).

## Quick Start

### Basic Usage
```bash
python code/generate_permutation_distribution.py
```

Or from the `code/` directory:
```bash
cd code
python generate_permutation_distribution.py
```

This will:
- Process all gene list sizes (50, 100, 150, ..., 1000)
- Run 1000 permutations per size
- Use all available CPU cores for parallel processing
- Resume from already completed permutations
- Save results to `results/permutation_results/`

### For c6i.8xlarge (32 vCPUs) - Recommended
```bash
python code/generate_permutation_distribution.py --n-jobs 32
```
This will complete in approximately **~3.5 hours** with maximum parallelization.

### Advanced Options

```bash
# Use specific number of CPU cores
python code/generate_permutation_distribution.py --n-jobs 8

# Process only specific gene list sizes
python code/generate_permutation_distribution.py --sizes 50 100 200

# Overwrite existing results (don't resume)
python code/generate_permutation_distribution.py --no-resume

# Change number of permutations per size
python code/generate_permutation_distribution.py --n-permutations 500
```

## Command Line Arguments

- `--n-permutations N`: Number of permutations per gene list size (default: 1000)
- `--n-jobs N`: Number of parallel jobs (default: number of CPU cores)
- `--resume`: Resume from already completed permutations (default: True)
- `--no-resume`: Do not resume, overwrite existing results
- `--sizes SIZE1 SIZE2 ...`: Process only specific gene list sizes (default: all sizes 50-1000)

## Output Structure

```
results/permutation_results/
├── config.json                    # Configuration used for this run
├── summary.json                   # Summary statistics
├── permutation_distribution.log   # Log file
├── size_50/
│   ├── permutation_0001.tsv
│   ├── permutation_0002.tsv
│   ├── ...
│   └── permutation_1000.tsv
├── size_100/
│   └── ...
└── size_1000/
    └── ...
```

## Output TSV Format

Each TSV file contains combined results from all 5 libraries with columns:
- `Library`: Library name (Reactome, KEGG, GO BP, GO MF, GO CC)
- `Iteration`: Iteration number
- `Term`: Term name
- `Description`: Term description
- `iteration overlapping genes`: Overlap size (e.g., "5/50")
- `iteration p-value`: P-value for this iteration
- `iteration -log(p-value)`: Negative log10 of p-value
- `Genes removed for next iteration`: Comma-separated list of genes
- Additional columns may include full list enrichment data

## Libraries Used

1. **Reactome**: `c2.cp.reactome.v2025.1.Hs.symbols.gmt`
2. **KEGG**: `c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt`
3. **GO Biological Process**: `c5.go.bp.v2025.1.Hs.symbols.gmt`
4. **GO Molecular Function**: `c5.go.mf.v2025.1.Hs.symbols.gmt`
5. **GO Cellular Component**: `c5.go.cc.v2025.1.Hs.symbols.gmt`

## Parameters

- **Background**: `all_genes.txt` from `data/backgrounds/`
- **P-value threshold**: 0.05 (raw p-value)
- **Min overlap**: 3
- **Min term size**: 10
- **Max term size**: 600
- **Max iterations**: 10
- **P-value method**: Fisher's Exact Test

## Performance

### Expected Runtime
- **Single-threaded**: ~110 hours total (5.5 hours per size × 20 sizes)
- **8 cores**: ~14 hours total
- **16 cores**: ~7 hours total

### Memory Requirements
- Libraries are loaded once and reused
- Background gene set is loaded once
- Each permutation processes one gene set at a time
- Estimated memory: 2-4 GB depending on library sizes

## Resumability

The script automatically detects and skips already completed permutations:
- Checks for existing TSV files in each size directory
- Extracts permutation index from filename
- Only processes missing permutations
- Use `--no-resume` to overwrite existing results

## Logging

Progress is logged to:
- **Console**: Real-time progress updates
- **File**: `permutation_distribution.log` (detailed logs)

## Error Handling

- Individual permutation failures are logged but don't stop the script
- Errors are recorded in the log file
- Summary statistics include success/failure counts
- Failed permutations can be re-run by resuming

## Example Workflow

```bash
# 1. Start the full run (will take many hours)
python code/generate_permutation_distribution.py --n-jobs 16

# 2. Check progress
tail -f permutation_distribution.log

# 3. If interrupted, resume (automatically skips completed)
python code/generate_permutation_distribution.py --n-jobs 16

# 4. Check summary
cat results/permutation_results/summary.json
```

## Troubleshooting

### Libraries Not Found
If you see warnings about missing library files:
- Check that all 5 library files exist in `data/libraries/`
- Verify filenames match exactly (case-sensitive)

### Background File Not Found
- Ensure `data/backgrounds/all_genes.txt` exists
- Check file permissions

### Out of Memory
- Reduce `--n-jobs` to use fewer parallel processes
- Process sizes one at a time using `--sizes`

### Slow Performance
- Increase `--n-jobs` if you have more CPU cores available
- Check disk I/O (writing 20,000 TSV files can be slow)
- Consider using faster storage (SSD)

## Notes

- Random seeds are based on permutation index for reproducibility
- Each permutation uses the same seed across runs
- Results are deterministic for the same permutation index
- Output files use TSV format (tab-separated values)

