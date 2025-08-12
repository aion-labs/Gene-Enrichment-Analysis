# Library Validation Module

## Overview

The `LibraryValidator` module ensures consistency between the gene_info database and MSigDB libraries by validating and updating gene symbols in GMT files. This addresses the issue where some genes in the libraries may not be present in the background gene set, leading to inconsistencies in enrichment analysis.

## Key Features

- **Gene Symbol Validation**: Uses the `GeneConverter` to validate gene symbols against the NCBI gene_info database
- **Automatic Updates**: Updates gene symbols to their current, validated versions
- **Safe Handling**: Keeps invalid genes unchanged rather than removing them
- **Backup Creation**: Automatically creates timestamped backups of original files
- **Batch Processing**: Can process entire directories of GMT files
- **Detailed Reporting**: Provides comprehensive statistics on validation results

## How It Works

### Gene Validation Process

1. **Direct Match**: Checks if the gene symbol exists directly in the gene_info database
2. **Synonyms**: Looks up the gene symbol in the synonyms database
3. **Gene History**: Maps old gene symbols to current symbols using gene history data
4. **Case Normalization**: Converts gene symbols to uppercase for consistency

### Update Strategy

- **Validated Genes**: Updated to their current, validated symbols
- **Invalid Genes**: Left unchanged (e.g., mitochondrial genes, Ensembl IDs)
- **Background Impact**: Invalid genes are automatically excluded during enrichment analysis since the background is library-specific

## Usage

### Command Line Interface

#### Update Libraries In Place (with backups)
```bash
python update_msigdb_libraries.py data/libraries/
```

#### Update Libraries and Save to New Directory
```bash
python update_msigdb_libraries.py data/libraries/ -o data/libraries_validated/
```

#### Update Without Creating Backups
```bash
python update_msigdb_libraries.py data/libraries/ --no-backup
```

#### Verbose Logging
```bash
python update_msigdb_libraries.py data/libraries/ -v
```

### Programmatic Usage

#### Single File Validation
```python
from code.library_validator import LibraryValidator

validator = LibraryValidator()
stats = validator.validate_gmt_file(
    "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt",
    create_backup=True
)
print(validator.get_validation_report())
```

#### Directory Validation
```python
from code.library_validator import LibraryValidator

validator = LibraryValidator()
stats = validator.validate_library_directory(
    "data/libraries/",
    output_dir="data/libraries_validated/",
    create_backup=True
)
print(validator.get_validation_report())
```

## Example Output

```
🧬 MSIGDB LIBRARY UPDATE
============================================================
Library directory: data/libraries/
Output directory: (overwrite original)
Create backups: True

🔧 Initializing library validator...
🔍 Validating and updating libraries...

============================================================
LIBRARY VALIDATION REPORT
============================================================
Files processed: 33
Terms processed: 15,847
Total genes: 975,544
Validated genes: 69 (0.1%)
Unchanged genes: 973,398 (99.9%)
Invalid genes: 77 (0.1%)
============================================================

📊 DETAILED ANALYSIS
----------------------------------------
Total genes processed: 975,544
Genes updated: 69 (0.01%)
Genes unchanged: 973,398 (99.89%)
Genes that could not be validated: 77 (0.01%)

⚠️  Note: 77 genes could not be validated.
   These genes will remain unchanged in the libraries.
   They will be ignored during enrichment analysis since the
   background is library-specific and will exclude them anyway.

✅ Library update completed successfully!
Original libraries have been updated in place.
Backup files have been created with timestamp suffixes.
```

## When to Use

### Recommended Usage Scenarios

1. **After Updating gene_info Database**: When you download a new version of the NCBI gene_info file
2. **After Updating MSigDB Libraries**: When you download new versions of MSigDB libraries
3. **Before Running Enrichment Analysis**: To ensure maximum consistency between libraries and background
4. **Quality Control**: To identify and fix gene symbol inconsistencies

### Frequency

- **Regular Updates**: Run after each gene_info or MSigDB update
- **Before Analysis**: Run before conducting new enrichment analyses
- **Quality Checks**: Run periodically to maintain consistency

## Benefits

### Statistical Accuracy
- **Improved Background Consistency**: Ensures library-specific backgrounds are more accurate
- **Better P-value Calculations**: Reduces discrepancies in hypergeometric test parameters
- **More Reliable Results**: Minimizes false positives/negatives due to gene symbol mismatches

### Data Quality
- **Standardized Gene Symbols**: All libraries use consistent, validated gene symbols
- **Reduced Missing Genes**: Minimizes genes that are in libraries but not in background
- **Better Documentation**: Clear tracking of which genes were updated and why

### Operational Efficiency
- **Automated Process**: No manual intervention required
- **Safe Updates**: Automatic backup creation prevents data loss
- **Comprehensive Reporting**: Detailed statistics on validation results

## File Structure

```
code/
├── library_validator.py          # Main validation module
├── gene_converter.py             # Gene symbol conversion utilities
└── ...

scripts/
├── update_msigdb_libraries.py    # Command-line interface
└── test_library_validator.py     # Test suite

data/
├── libraries/                    # MSigDB GMT files
│   ├── c2.cp.reactome.v2025.1.Hs.symbols.gmt
│   ├── c5.go.bp.v2025.1.Hs.symbols.gmt
│   └── ...
└── backgrounds/
    └── all_genes.txt            # Background gene set
```

## Troubleshooting

### Common Issues

1. **No GMT Files Found**: Ensure the library directory contains `.gmt` files
2. **Permission Errors**: Check file permissions for read/write access
3. **Memory Issues**: For very large libraries, process files individually
4. **Backup Failures**: Ensure sufficient disk space for backup files

### Validation Failures

- **Mitochondrial Genes**: MT-* genes often cannot be validated (expected)
- **Ensembl IDs**: ENSG* identifiers are not gene symbols (expected)
- **Novel Genes**: Recently discovered genes may not be in gene_info yet
- **Species Mismatch**: Ensure gene_info file is for the correct species

## Integration with Workflow

### Pre-Analysis Setup
```bash
# 1. Update gene_info database
# 2. Update MSigDB libraries
# 3. Validate libraries for consistency
python update_msigdb_libraries.py data/libraries/

# 4. Run enrichment analysis
python -m streamlit run code/streamlit_app.py
```

### Automated Pipeline
```bash
#!/bin/bash
# Example automated workflow

# Update gene_info
wget -O Homo_sapiens.gene_info.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
gunzip Homo_sapiens.gene_info.gz

# Update MSigDB libraries
# ... download new MSigDB files ...

# Validate libraries
python update_msigdb_libraries.py data/libraries/ -o data/libraries_validated/

# Update alias.json if needed
# ... update library metadata ...

echo "Library validation completed successfully!"
```

## Conclusion

The LibraryValidator module provides a robust solution for maintaining consistency between gene annotation databases and MSigDB libraries. By automatically validating and updating gene symbols, it ensures more accurate enrichment analysis results while preserving data integrity through safe update practices.
