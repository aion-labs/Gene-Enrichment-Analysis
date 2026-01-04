# Reorganization Summary

This document summarizes the reorganization completed for the release preparation.

## New Folder Structure

### `tools/` - User-Accessible Utilities
Contains tools for:
- **Generating permutations**: `generate_permutation_distribution.py`, `generate_permutation_statistics_report.py`
- **Managing libraries**: `update_msigdb.py`, `update_msigdb_simple.py`, `update_msigdb_libraries.py`, `create_stringdb_library.py`
- **Preparing benchmarking data**: `extract_cluster_statistics_from_permutations.py`, `prepare_parquet_cluster_stats.py`, `unzip_and_merge_permutations.py`
- **Validating libraries**: `library_validator.py`
- **Documentation**: `README.md` with usage instructions

### `development/` - Development and Testing
- `tests/` - All test files (`test_*.py`)
- `scripts/` - Development/debug scripts
- `test_results/` - Test output directories
- `data/` - Development data files (logs, statistics)

### `code/` - Core Application Only
Contains only the essential application code:
- CLI (`cli.py`)
- Streamlit app (`streamlit_app.py`)
- Core enrichment logic (`enrichment.py`, `iter_enrichment.py`)
- Supporting modules (`background_gene_set.py`, `gene_set.py`, `gene_set_library.py`, `gene_converter.py`)
- Benchmarking core (`network_connectivity_benchmark.py`, `parallel_null_distribution.py`)
- UI components (`ui/` folder)
- Static assets (`static/` folder)

### `archive/` - Historical Artifacts
- `documentation/` - Development documentation
- `results/` - Old test runs and outputs
- `permutations/` - Old permutation archives
- `data/` - Old library backups
- `permutation_analysis/` - Existing archive (preserved)

## Files Moved

### From `code/` to `tools/`:
- `generate_permutation_distribution.py`
- `generate_permutation_statistics_report.py`
- `library_validator.py`
- `run_permutations.sh` (if existed)

### From root to `tools/`:
- `update_msigdb.py`
- `update_msigdb_simple.py`
- `update_msigdb_libraries.py`
- `create_stringdb_library.py`

### From `permutations/` to `tools/`:
- `extract_cluster_statistics_from_permutations.py`
- `prepare_parquet_cluster_stats.py`
- `unzip_and_merge_permutations.py`

### From root to `development/tests/`:
- All `test_*.py` files
- `test_all_libraries/` folder

### From root to `development/scripts/`:
- Analysis scripts (`analyze_*.py`, `benchmark_*.py`, etc.)
- Debug scripts (`debug_*.py`)
- Example scripts (`example_*.py`)
- Extraction/merge scripts
- Development shell scripts

### From root to `development/test_results/`:
- `cli_test_results/`
- `cli_test_results_benchmark/`
- `cli_results/`

### From root to `archive/documentation/`:
- Development documentation files (see plan for full list)

### From root to `archive/results/`:
- `results/` folder (old test runs)

### From `permutations/` to `archive/permutations/`:
- Old permutation result archives (`results-*.zip`)

## Files Removed

### From `code/`:
- `sandbox.py` (debug file)
- `streamlit_app_test.py` (test version)

## Import Path Updates

All moved files have been updated to use correct import paths:
- Tools in `tools/` now import from `code/` using `sys.path.insert(0, str(PROJECT_ROOT / "code"))`
- `PROJECT_ROOT` paths updated to account for new folder locations
- All imports verified to work correctly

## User-Facing Documentation (Kept in Root)

The following documentation files remain in the root as they are user-facing:
- `README.md` - Main documentation
- `CLI_README.md` - CLI documentation
- `ENTREZ_ID_SUPPORT.md` - User guide for Entrez ID support
- `LIBRARY_VALIDATION_README.md` - User guide for library validation
- `LINUX_INSTALLATION.md` - Installation guide
- `MSIGDB_UPDATE_README.md` - Guide for updating MSigDB
- `REPRODUCING.md` - Reproducibility guide
- `STREAMLIT_CLOUD_DEPLOYMENT.md` - Deployment guide
- `CHANGELOG.md` - Release notes
- `LICENSE` - License file

## Data Files

- All data files remain in `data/` folder
- `Homo_sapiens.gene_info` exists in both root and `data/recent_release/` (files differ, both kept)
- Benchmarking data remains in `permutations/` folder

## Verification

The reorganization maintains:
- ✅ All core application functionality
- ✅ All tools accessible and functional
- ✅ Clear separation between user tools and development code
- ✅ Easy access to library management and permutation generation tools
- ✅ Clean code folder with only essential application code

## Next Steps

1. Test CLI: `python code/cli.py --help`
2. Test Streamlit: `streamlit run code/streamlit_app.py`
3. Test tools: `python tools/library_validator.py --help`
4. Verify benchmarking works with existing permutation data
5. Review user-facing documentation files and move to `tools/` if they're tool-specific

