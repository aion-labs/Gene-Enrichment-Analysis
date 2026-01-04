# Installation Tests

This directory contains tests to verify that the installation is complete and functional.

## Test Scripts

### `test_installation.py`
**Comprehensive installation verification test**

Runs all installation checks:
- Python version verification (>= 3.12)
- Dependency checking (all required packages)
- File structure verification
- CLI functionality test
- Streamlit functionality test

**Usage:**
```bash
python test/test_installation.py
```

### `test_cli.py`
**CLI functionality test**

Tests the command-line interface:
- Regular enrichment mode
- Iterative enrichment mode
- Output file generation
- Result validation

**Usage:**
```bash
python test/test_cli.py
```

### `test_streamlit.py`
**Streamlit app installation test**

Verifies that:
- Streamlit and dependencies can be imported
- The streamlit_app module loads correctly
- Required data files are accessible
- Static files (logos) are available

**Usage:**
```bash
python test/test_streamlit.py
```

### `test_paths.py`
**Path resolution test**

Verifies that file paths resolve correctly:
- Gene info file paths
- Data directory paths
- Code directory paths

**Usage:**
```bash
python test/test_paths.py
```

## Running All Tests

To run all installation tests:

```bash
# From project root
python test/test_installation.py
```

This will run all checks and provide a comprehensive report.

## Expected Results

After a successful installation, all tests should pass:
- ✅ Python version check
- ✅ All dependencies installed
- ✅ File structure correct
- ✅ CLI works correctly
- ✅ Streamlit app can be imported

## Troubleshooting

If tests fail:

1. **Dependencies missing**: Install with `pip install -r requirements.txt`
2. **File structure issues**: Ensure you're in the project root directory
3. **Path issues**: Check that data files are in the correct locations
4. **Import errors**: Verify Python path includes the `code/` directory

## Notes

- These tests are designed to be run from the project root directory
- Some tests may take a few minutes (especially CLI tests that run actual analyses)
- Test output files are written to the `sandbox/` directory

