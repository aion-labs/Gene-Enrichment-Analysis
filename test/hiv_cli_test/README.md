# HIV Dependency Factors CLI Test

This directory contains a test script that runs the HIV dependency factors gene list through the CLI with all 11 libraries that have permutation data available.

## Test Configuration

- **Gene List**: `HIV_dependency_factors.symbols.txt` (272 genes)
- **Background**: `all_genes.txt`
- **Libraries**: All 11 libraries with permutation data:
  1. Reactome (c2.cp.reactome.v2025.1.Hs.symbols.gmt)
  2. KEGG (c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt)
  3. GO BP (c5.go.bp.v2025.1.Hs.symbols.gmt)
  4. GO MF (c5.go.mf.v2025.1.Hs.symbols.gmt)
  5. GO CC (c5.go.cc.v2025.1.Hs.symbols.gmt)
  6. BioCarta (c2.cp.biocarta.v2025.1.Hs.symbols.gmt)
  7. Canonical pathways (c2.cp.v2025.1.Hs.symbols.gmt)
  8. KEGG Medicus (c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt)
  9. Pathway Interaction Database (c2.cp.pid.v2025.1.Hs.symbols.gmt)
  10. WikiPathways (c2.cp.wikipathways.v2025.1.Hs.symbols.gmt)
  11. Hallmark (h.all.v2025.1.Hs.symbols.gmt)

- **P-value threshold**: 0.01
- **Mode**: Iterative (with benchmarking enabled)
- **Other parameters**: Default values

## Usage

From the project root:

```bash
# Activate virtual environment first
source .venv/bin/activate

# Run the test
python test/hiv_cli_test/run_hiv_test.py
```

Or directly:

```bash
python test/hiv_cli_test/run_hiv_test.py
```

## Expected Output

The test will:
1. Check prerequisites (gene list, background, libraries, CLI)
2. Run the CLI analysis with all 11 libraries
3. Generate results in `test/hiv_cli_test/results/`
4. Verify that expected output files were created

## Output Files

Results will be saved to `test/hiv_cli_test/results/HIV_dependency_factors_<timestamp>/`:

- `combined_iterative_results.tsv` - Combined results from all libraries
- `*_iterative_results.tsv` - Individual library results (one per library)
- `iterative_enrichment_snapshot.json` - Metadata and summary
- `combined_network.dot` - Network visualization file
- `statistical_benchmarks.json` - Benchmark statistics (if benchmarking enabled)
- `statistical_benchmarks_table.tsv` - Benchmark table (if benchmarking enabled)

## Notes

- This test may take several minutes to complete (depending on system resources)
- Benchmarking requires permutation data to be available in `permutations/permutation_cluster_statistics_parquet/`
- The test uses iterative mode to enable statistical benchmarking
- All 11 libraries are used to match the permutation data structure

