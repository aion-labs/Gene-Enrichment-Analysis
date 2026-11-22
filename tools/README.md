# Tools for Gene Enrichment Analysis

This directory contains utilities for managing gene set libraries, generating permutation data for benchmarking, and preparing data files.

## Table of Contents

- [Generating Permutations](#generating-permutations)
- [Managing Libraries](#managing-libraries)
- [Preparing Benchmarking Data](#preparing-benchmarking-data)
- [Validating Libraries](#validating-libraries)

## Generating Permutations

### `generate_permutation_distribution.py`

Generates random gene lists and runs iterative enrichment analysis to create a null distribution for statistical benchmarking.

**Usage:**
```bash
python tools/generate_permutation_distribution.py [options]
```

**Options:**
- `--n-permutations N`: Number of permutations to generate (default: 1000)
- `--n-jobs N`: Number of parallel jobs (default: number of CPUs)
- `--resume`: Resume from previous run
- `--gene-list-sizes SIZES`: Comma-separated list of gene list sizes (default: 50,100,150,...,1000)

**Example:**
```bash
# Generate 1000 permutations using 32 parallel jobs
python tools/generate_permutation_distribution.py --n-permutations 1000 --n-jobs 32
```

**Output:**
- Permutation results saved to `results/permutation_results/`
- Each permutation generates a TSV file with enrichment results

### `generate_permutation_statistics_report.py`

Generates statistical reports from permutation data across different gene list sizes and p-value thresholds.

**Usage:**
```bash
python tools/generate_permutation_statistics_report.py [options]
```

**Output:**
- Summary report (text file)
- Detailed CSV with cluster statistics

## Managing Libraries

Use **`code/libraries.py`** (see the main README, "Managing gene-set collections") to list, activate,
deactivate, add and remove collections, and to update MSigDB and the NCBI gene data. It harmonises gene
symbols with `data/recent_release` exactly as the bundled MSigDB 2025.1 files were prepared (re-running it on
the official 2025.1 zip reproduces them byte for byte), keeps collection names, colours and active flags, and
moves replaced files to `data/libraries/previous/` instead of deleting them.

The older scripts below are kept for reference only. Do not use them for analyses reported with this
version: `update_msigdb.py` copies only C2 and C5 files (not Hallmark), gives collections different names
from `alias.json`, adds new files inactive alongside old versions, and **rewrites `data/backgrounds/all_genes.txt`**
from the C2/C5 genes, which changes every result; it also takes the MSigDB password on the command line.
`update_msigdb_libraries.py` / `library_validator.py` perform the symbol harmonisation now built into
`code/libraries.py` (and require Python 3.12+).

### `create_stringdb_library.py`

Creates a gene set library from STRING-DB protein interactions.

**Usage:**
```bash
python tools/create_stringdb_library.py
```

**What it does:**
1. Downloads STRING-DB human protein interaction file
2. Filters interactions by score >= 0.9
3. Maps Ensembl IDs to gene symbols
4. Creates GMT format library
5. Adds to the enrichment system

**Requirements:**
- `Homo_sapiens.gene_info` file in `data/recent_release/`
- Internet connection for downloading STRING-DB data

## Preparing Benchmarking Data

### `extract_cluster_statistics_from_permutations.py`

Extracts cluster-level statistics from permutation results.

**Usage:**
```bash
python tools/extract_cluster_statistics_from_permutations.py [options]
```

**Output:**
- TSV file with cluster statistics for each permutation
- Columns include: cluster size, number of genes/terms, density, libraries, etc.

### `prepare_parquet_cluster_stats.py`

Prepares cluster statistics in Parquet format, segmented by gene list size.

**Usage:**
```bash
python tools/prepare_parquet_cluster_stats.py
```

**Output:**
- Parquet files in `permutations/permutation_cluster_statistics_parquet/`
- One file per gene list size: `cluster_stats_size_*.parquet`

**Note:** These Parquet files are used by the benchmarking system for fast loading.

### `unzip_and_merge_permutations.py`

Processes permutation result archives and merges them into a single file.

**Usage:**
```bash
python tools/unzip_and_merge_permutations.py [options]
```

## Validating Libraries

### `library_validator.py`

Validates gene symbols in GMT files against the NCBI gene_info database and updates them to use consistent symbols.

**Usage:**
```bash
# Validate a single file
python tools/library_validator.py path/to/library.gmt -o path/to/output.gmt

# Validate all files in a directory
python tools/library_validator.py data/libraries/ -o data/libraries_validated/

# Validate without creating backups
python tools/library_validator.py data/libraries/ --no-backup
```

**Options:**
- `-o, --output`: Output file or directory
- `--no-backup`: Don't create backup files
- `--pattern PATTERN`: File pattern for directory processing (default: `*.gmt`)
- `--gene-info PATH`: Path to Homo_sapiens.gene_info file
- `--gene-history PATH`: Path to gene_history file
- `-v, --verbose`: Verbose logging

**What it does:**
- Validates gene symbols against NCBI gene_info database
- Updates old gene symbols to current official symbols
- Reports validation statistics
- Creates backups of original files (unless `--no-backup` is used)

## Common Workflows

### Adding a New Library

1. **Download or create your GMT file** in the standard GMT format
2. **Validate the library:**
   ```bash
   python tools/library_validator.py your_library.gmt -o data/libraries/your_library.gmt
   ```
3. **Add to alias.json** in `data/libraries/` to make it available in the application

### Generating Benchmarking Data

1. **Generate permutations:**
   ```bash
   python tools/generate_permutation_distribution.py --n-permutations 1000 --n-jobs 32
   ```
2. **Extract cluster statistics:**
   ```bash
   python tools/extract_cluster_statistics_from_permutations.py
   ```
3. **Prepare Parquet files:**
   ```bash
   python tools/prepare_parquet_cluster_stats.py
   ```
4. **Generate statistics report:**
   ```bash
   python tools/generate_permutation_statistics_report.py
   ```

### Updating MSigDB Libraries

1. **Download latest MSigDB release** from https://www.gsea-msigdb.org/gsea/msigdb
2. **Update libraries:**
   ```bash
   python tools/update_msigdb.py --zip-file path/to/msigdb_v2025.1.Hs_files_to_download_locally.zip
   ```
3. **Validate updated libraries:**
   ```bash
   python tools/update_msigdb_libraries.py data/libraries/
   ```

## Requirements

All tools require:
- Python 3.12+
- Dependencies from `requirements.txt`
- Access to `data/` directory with libraries and backgrounds
- For some tools: `Homo_sapiens.gene_info` in `data/recent_release/`

## Notes

- All tools are designed to be run from the project root directory
- Tools that modify libraries create backups by default
- Check tool help messages for detailed options: `python tools/<script>.py --help`
- Some tools require significant computational resources (especially permutation generation)

