# Release Preparation Plan

## Overview
This document outlines the plan to organize the codebase for release with a clear structure:
- **`code/`** - Core application code (CLI, Streamlit, enrichment logic)
- **`tools/`** - Utilities for generating permutations and managing libraries (user-accessible)
- **`development/`** - Development code and tests
- **`archive/`** - Historical documentation, old results, and other artifacts

## Summary of Changes

### New Folders to Create
1. **`tools/`** - Contains utilities for:
   - Generating permutation data for benchmarking
   - Adding/updating gene set libraries
   - Preparing and processing data files
   - Includes a `README.md` with usage instructions

2. **`development/`** - Contains:
   - Test files (`tests/`)
   - Development/debug scripts (`scripts/`)
   - Test results (`test_results/`)
   - Development data files (`data/`)

### Files to Move

**From `code/` to `tools/`:**
- `generate_permutation_distribution.py`
- `generate_permutation_statistics_report.py`
- `library_validator.py`

**From root to `tools/`:**
- `update_msigdb.py`, `update_msigdb_simple.py`, `update_msigdb_libraries.py`
- `create_stringdb_library.py`

**From `permutations/` to `tools/`:**
- `extract_cluster_statistics_from_permutations.py`
- `prepare_parquet_cluster_stats.py`
- `unzip_and_merge_permutations.py`

**From root to `development/tests/`:**
- All `test_*.py` files
- `test_all_libraries/` folder

**From root to `development/scripts/`:**
- Analysis, debug, and example scripts
- Development shell scripts

**From root to `archive/documentation/`:**
- Most `.md` files (except `README.md`, `CLI_README.md`, and user-facing docs)

**From root to `archive/results/`:**
- `results/` folder (old test runs)

**From root to `archive/`:**
- Old permutation archives
- Old library backups

## Essential Components (Keep in Root/Code/Data)

### 1. Core Application Code (`code/` folder)
**Keep all files:**
- `cli.py` - Command-line interface
- `streamlit_app.py` - Web interface
- `enrichment.py` - Core enrichment logic
- `iter_enrichment.py` - Iterative enrichment logic
- `background_gene_set.py` - Background gene set handling
- `gene_set.py` - Gene set handling
- `gene_set_library.py` - Library management
- `gene_converter.py` - Gene ID conversion
- `network_connectivity_benchmark.py` - Benchmarking core
- `parallel_null_distribution.py` - Null distribution computation
- `ui/` folder - All UI components for Streamlit
- `static/` folder - Logo images
- `entry.sh` - Docker entry point
- `run` - Execution script (if needed)

**Remove from code/:**
- `sandbox.py` - Debug file (already in .gitignore)
- `streamlit_app_test.py` - Test version of streamlit app
- `generate_permutation_distribution.py` - Move to `tools/`
- `generate_permutation_statistics_report.py` - Move to `tools/`
- `library_validator.py` - Move to `tools/`
- `run_permutations.sh` - Move to `tools/` (if needed)

### 2. Data Files (`data/` folder)
**Keep all:**
- `backgrounds/` - All background gene lists
- `libraries/` - All GMT library files and alias.json
- `gene_lists/` - Example gene lists
- `recent_release/` - NCBI gene info files (needed for Entrez ID support)

**Note:** `data/libraries/backup/` contains many backup folders - consider archiving old backups

### 3. Benchmarking Data (`permutations/` folder)
**Keep:**
- `permutation_cluster_statistics_parquet/` - Parquet files needed for benchmarking
- `merged_permutation_results.tsv` - Merged permutation data
- `permutation_cluster_statistics.tsv` - Cluster statistics

**Move to `tools/`:**
- `extract_cluster_statistics_from_permutations.py` - Extract statistics from permutations
- `prepare_parquet_cluster_stats.py` - Prepare parquet files for benchmarking
- `unzip_and_merge_permutations.py` - Process permutation archives

**Archive:**
- `results-*.zip` - Old permutation result archives

### 4. Root Level Files
**Keep:**
- `README.md` - Main documentation
- `CLI_README.md` - CLI documentation
- `LICENSE` - License file
- `requirements.txt` - Python dependencies
- `pyproject.toml` - Project configuration
- `run_example.sh` - Example usage script
- `install_linux.sh` - Installation script
- `packages.txt` - System dependencies (if needed)

**Move to archive:**
- All other `.md` files (development documentation)
- `Homo_sapiens.gene_info` - Should be in `data/recent_release/` (check if duplicate)

## Tools Folder (`tools/`)

Create a new `tools/` folder for utilities that help users:
- Generate permutation data for benchmarking
- Add or update gene set libraries
- Prepare and process data files

### Tools for Permutation Generation
Move from `code/` to `tools/`:
- `generate_permutation_distribution.py` - Generate permutation data
- `generate_permutation_statistics_report.py` - Generate statistics reports

Move from `permutations/` to `tools/`:
- `extract_cluster_statistics_from_permutations.py` - Extract statistics
- `prepare_parquet_cluster_stats.py` - Prepare parquet files
- `unzip_and_merge_permutations.py` - Process permutation archives

### Tools for Library Management
Move from root to `tools/`:
- `update_msigdb.py` - Update MSigDB libraries
- `update_msigdb_simple.py` - Simple MSigDB update script
- `update_msigdb_libraries.py` - Library update utility
- `create_stringdb_library.py` - Create STRING-DB library
- `library_validator.py` - Validate library files (from `code/`)

### Tools Documentation
**Create `tools/README.md`** with clear instructions on:
- **Generating Permutations**: How to use `generate_permutation_distribution.py` to create permutation data for benchmarking
- **Adding Libraries**: How to use `create_stringdb_library.py` and other library creation tools
- **Updating Libraries**: How to use `update_msigdb.py` and related scripts to update MSigDB libraries
- **Preparing Benchmarking Data**: How to use extraction and preparation scripts to process permutation results
- **Validating Libraries**: How to use `library_validator.py` to validate gene set libraries

This README should make it easy for users to:
1. Generate their own permutation data for custom benchmarking
2. Add custom gene set libraries
3. Update existing libraries with new versions

## Development Folder (`development/`)

Create a new `development/` folder for:
- Test files
- Development/debug scripts
- Test results

### 1. Test Files
Move to `development/tests/`:
- `test_*.py` (all test files from root)
- `test_all_libraries/` folder

### 2. Development Scripts
Move to `development/scripts/`:
- `analyze_library_coverage.py`
- `benchmark_connectivity_computation.py`
- `benchmark_hiv_connectivity.py`
- `check_reactome_genes.py`
- `combine_permutation_data.py`
- `compress_old_results.py`
- `debug_network.py`
- `example_ai_analysis_formats.py`
- `example_network_benchmark.py`
- `extract_cluster_stats_from_combined_data.py`
- `extract_cluster_stats_from_firstrun.py`
- `find_missing_genes.py`
- `generate_average_cluster_networks.py`
- `generate_cluster_statistical_report.py`
- `library_statistics.py`
- `merge_all_permutation_files.py`
- `merge_combined_permutation_data.py`
- `rebuild_parquet_benchmarking.sh`
- `launch_gene_enrichment.sh`
- `Gene_Enrichment.command`

### 3. Test Results
Move to `development/test_results/`:
- `cli_test_results/`
- `cli_test_results_benchmark/`
- `cli_results/`

### 4. Development Data Files
Move to `development/data/`:
- `library_statistics_output.txt`
- `library_statistics_summary.csv`
- `benchmarking.log`
- `extract_cluster_stats.log`

## Archive Folder (`archive/`)

Move historical and non-essential items to `archive/`:

### 1. Development Documentation
Move to `archive/documentation/`:
- `AI_ANALYSIS_FORMATS.md`
- `analyze_job_log.md`
- `AWS_INSTANCE_RECOMMENDATIONS.md`
- `BENCHMARK_WORKFLOW_EXPLANATION.md`
- `BENCHMARKING_METHODOLOGY_CONFIRMED.md`
- `BENCHMARKING_WITH_RANDOM_PERMUTATIONS_EXPLAINED.md`
- `CHANGELOG.md` (or keep if needed for release notes)
- `CLUSTER_BASED_BENCHMARKING.md`
- `CLUSTER_BASED_CONNECTIVITY_IMPLEMENTATION.md`
- `CLUSTER_BASED_CONNECTIVITY_PROPOSAL.md`
- `CLUSTER_BENCHMARKING_APPROACH.md`
- `CLUSTER_BENCHMARKING_FIXES.md`
- `CONNECTIVITY_SIZE_HANDLING.md`
- `DENSITY_CALCULATION_EXPLAINED.md`
- `DENSITY_NORMALIZATION_PROPOSAL.md`
- `ENTREZ_ID_SUPPORT.md` (or keep if user-facing)
- `HYPERGEOMETRIC_VERIFICATION.md`
- `LIBRARY_VALIDATION_README.md` (or keep if user-facing)
- `LINUX_INSTALLATION.md` (or keep if user-facing)
- `LOGGING_IMPROVEMENTS.md`
- `MSIGDB_UPDATE_README.md` (or keep if user-facing)
- `NETWORK_BENCHMARK_SUMMARY.md`
- `NETWORK_CONNECTIVITY_BENCHMARKING.md`
- `NULL_DISTRIBUTION_METHODOLOGY_VERIFICATION.md`
- `NULL_DISTRIBUTION_METRICS_FIX.md`
- `PARALLEL_NULL_DISTRIBUTION_IMPLEMENTATION.md`
- `PERMUTATION_STATISTICS_REPORT_USAGE.md`
- `REGENERATE_PERMUTATIONS.md`
- `REPRODUCING.md` (or keep if user-facing)
- `STATISTICAL_ANALYSIS_EXPLAINED.md`
- `STREAMLIT_CLOUD_DEPLOYMENT.md` (or keep if user-facing)

### 2. Results Folder
Move to `archive/results/`:
- `results/` folder (contains old test runs and development outputs)

### 3. Old Permutation Archives
Move to `archive/permutations/`:
- `permutations/results-*.zip` (old permutation result archives)

### 4. Library Backups
Move to `archive/data/library_backups/`:
- Old backup folders from `data/libraries/backup/` (keep only most recent if needed)

## Data Organization

### Homo_sapiens.gene_info Location
- Check if `Homo_sapiens.gene_info` exists in root
- If yes, verify it's the same as `data/recent_release/Homo_sapiens.gene_info`
- Remove root copy if duplicate (code expects it in `data/recent_release/`)

### Library Backups
- `data/libraries/backup/` contains many timestamped backup folders
- Consider archiving old backups to `archive/data/library_backups/`
- Keep only the most recent backup if needed

## Final Folder Structure

```
project_root/
├── code/                          # Core application
│   ├── cli.py
│   ├── streamlit_app.py
│   ├── enrichment.py
│   ├── iter_enrichment.py
│   ├── ui/
│   └── ...
├── tools/                         # Utilities for users
│   ├── README.md                  # How to use tools
│   ├── generate_permutation_distribution.py
│   ├── generate_permutation_statistics_report.py
│   ├── update_msigdb.py
│   ├── create_stringdb_library.py
│   ├── library_validator.py
│   └── ...
├── development/                   # Development and testing
│   ├── tests/                     # Test files
│   ├── scripts/                   # Dev/debug scripts
│   ├── test_results/              # Test outputs
│   └── data/                      # Dev data files
├── data/                          # Application data
│   ├── backgrounds/
│   ├── libraries/
│   ├── gene_lists/
│   └── recent_release/
├── permutations/                  # Benchmarking data
│   ├── permutation_cluster_statistics_parquet/
│   ├── merged_permutation_results.tsv
│   └── permutation_cluster_statistics.tsv
├── archive/                       # Historical artifacts
│   ├── permutation_analysis/     # (already exists)
│   ├── documentation/             # Old docs
│   ├── results/                   # Old results
│   ├── permutations/              # Old permutation archives
│   └── data/                      # Old data backups
├── README.md
├── CLI_README.md
├── LICENSE
├── requirements.txt
├── pyproject.toml
├── run_example.sh
└── install_linux.sh
```

## Verification Checklist

After reorganization, verify:
- [ ] CLI runs: `python code/cli.py --help`
- [ ] Streamlit runs: `streamlit run code/streamlit_app.py`
- [ ] Benchmarking works (if permutation data exists)
- [ ] All imports resolve correctly in `code/` (no broken imports)
- [ ] Tools can be run from `tools/` directory (test imports work)
- [ ] Example script runs: `./run_example.sh`
- [ ] Data files are accessible from code
- [ ] `tools/README.md` exists and is clear
- [ ] No broken paths in documentation
- [ ] `Homo_sapiens.gene_info` is in correct location (`data/recent_release/`)

## Documentation Files

### Keep in Root (User-Facing)
- `README.md` - Main documentation
- `CLI_README.md` - CLI documentation

### Move to `tools/` (If User-Facing)
- `MSIGDB_UPDATE_README.md` - If users need to update libraries
- `LIBRARY_VALIDATION_README.md` - If users need to validate libraries

### Move to `archive/documentation/` (Development Docs)
- `AI_ANALYSIS_FORMATS.md`
- `analyze_job_log.md`
- `AWS_INSTANCE_RECOMMENDATIONS.md`
- `BENCHMARK_WORKFLOW_EXPLANATION.md`
- `BENCHMARKING_METHODOLOGY_CONFIRMED.md`
- `BENCHMARKING_WITH_RANDOM_PERMUTATIONS_EXPLAINED.md`
- `CHANGELOG.md` (or keep if needed for release notes)
- `CLUSTER_BASED_BENCHMARKING.md`
- `CLUSTER_BASED_CONNECTIVITY_IMPLEMENTATION.md`
- `CLUSTER_BASED_CONNECTIVITY_PROPOSAL.md`
- `CLUSTER_BENCHMARKING_APPROACH.md`
- `CLUSTER_BENCHMARKING_FIXES.md`
- `CONNECTIVITY_SIZE_HANDLING.md`
- `DENSITY_CALCULATION_EXPLAINED.md`
- `DENSITY_NORMALIZATION_PROPOSAL.md`
- `ENTREZ_ID_SUPPORT.md` (review - might be user-facing)
- `HYPERGEOMETRIC_VERIFICATION.md`
- `LINUX_INSTALLATION.md` (review - might be user-facing)
- `LOGGING_IMPROVEMENTS.md`
- `NETWORK_BENCHMARK_SUMMARY.md`
- `NETWORK_CONNECTIVITY_BENCHMARKING.md`
- `NULL_DISTRIBUTION_METHODOLOGY_VERIFICATION.md`
- `NULL_DISTRIBUTION_METRICS_FIX.md`
- `PARALLEL_NULL_DISTRIBUTION_IMPLEMENTATION.md`
- `PERMUTATION_STATISTICS_REPORT_USAGE.md`
- `REGENERATE_PERMUTATIONS.md`
- `REPRODUCING.md` (review - might be user-facing)
- `STATISTICAL_ANALYSIS_EXPLAINED.md`
- `STREAMLIT_CLOUD_DEPLOYMENT.md` (review - might be user-facing)

**Decision:** Review user-facing docs and keep in root or move to `tools/` if they help users add libraries or generate permutations.

