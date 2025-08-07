# Entrez ID Support

This document describes the new Entrez ID support feature added to the gene set enrichment analysis application.

## Overview

The application now supports input of gene identifiers in two distinct formats:
- **Gene symbols** (e.g., TP53, BRCA1)
- **Entrez IDs** (e.g., 7157, 672)

Users must explicitly specify the input format before entering their gene list.

**Note**: Both input gene lists and background gene lists support Entrez ID conversion.

## How It Works

### Gene Database
The application uses the NCBI gene information database (`Homo_sapiens.gene_info`) which contains:
- **193,593 Entrez ID to symbol mappings**
- **69,767 synonym to symbol mappings**
- Only human genes (taxonomy ID: 9606)

### Conversion Process
1. **Format selection**: User explicitly selects input format (Gene Symbols or Entrez IDs)
2. **Input parsing**: The application parses user input line by line
3. **Format-specific processing**:
   - **Entrez IDs**: Converted to official gene symbols
   - **Gene Symbols**: Validated against the database
4. **Validation**: Unrecognized identifiers are flagged for user review
5. **Analysis**: Only successfully processed gene symbols are used in enrichment analysis

## Usage

### Input Format Selection
Users must first select their input format:

1. **Gene Symbols**: Enter official gene symbols
   ```
   TP53
   BRCA1
   BRCA2
   ```

2. **Entrez IDs**: Enter numeric Entrez IDs
   ```
   7157
   672
   675
   ```

### Example Input
The "Input an example" button provides format-specific examples:

**Gene Symbols Example:**
```
ABHD11
ABHD14A
ABHD3
ACAA1A
ACBD4
ACO1
ADH5
ADHFE1
ADK
AFAP1L1
```

**Entrez IDs Example:**
```
7157
672
675
7159
7160
7158
7161
7162
7163
7164
```

### Background Gene List Support
Background gene lists also support both formats:

1. **Upload Format Selection**: Choose between "Gene Symbols" or "Entrez IDs" for background files
2. **Automatic Conversion**: Entrez IDs in background files are converted to symbols
3. **Validation**: Gene symbols in background files are validated against the database
4. **Feedback**: Success/error messages show the number of genes loaded

### Validation Results
The application displays a validation summary showing:
- ✅ **Valid genes** (successfully processed)
- ⚠️ **Invalid Entrez IDs** (not found in database) - when using Entrez ID format
- ⚠️ **Invalid symbols** (not found in database) - when using Symbol format
- 📊 **Database statistics** (total mappings available)

## Technical Implementation

### Files Added/Modified
- `code/gene_converter.py` - New gene conversion module
- `code/ui/helpers.py` - Updated with conversion functions
- `code/streamlit_app.py` - Integrated conversion into main app
- `code/background_gene_set.py` - Updated to support Entrez ID conversion
- `Homo_sapiens.gene_info` - NCBI gene database file

### Key Classes
- `GeneConverter`: Main conversion class
  - `convert_input()`: Converts mixed input to symbols
  - `get_symbol()`: Gets symbol for Entrez ID
  - `get_entrez()`: Gets Entrez ID for symbol
  - `get_stats()`: Returns database statistics

### Database Structure
The NCBI gene info file contains tab-separated fields:
- Column 1: Taxonomy ID (9606 for human)
- Column 2: Entrez ID
- Column 3: Official gene symbol
- Column 5: Synonyms (pipe-separated)

## Benefits

1. **Flexibility**: Users can input genes in their preferred format
2. **Compatibility**: Works with data from various sources
3. **Validation**: Clear feedback on unrecognized identifiers
4. **Comprehensive**: Includes synonyms for better matching
5. **Performance**: Fast lookup using dictionary-based mapping

## Testing

Run the test script to verify functionality:
```bash
python3 test_gene_converter.py
```

This will test:
- Database loading
- Mixed format conversion
- Individual ID lookups
- Error handling

## Data Source

The gene database is sourced from NCBI:
- **URL**: https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
- **Format**: Tab-separated text file
- **Update frequency**: Regular updates from NCBI
- **License**: Public domain (NCBI data)

## Future Enhancements

Potential improvements could include:
- Support for other species
- Additional identifier types (Ensembl IDs, RefSeq, etc.)
- Automatic database updates
- Custom synonym mappings
- Batch file processing 