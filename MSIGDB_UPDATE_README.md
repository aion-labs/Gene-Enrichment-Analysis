# MSigDB Update Scripts

This directory contains scripts to update the MSigDB gene set libraries used by the gene enrichment analysis application.

## Overview

The scripts can update MSigDB gene set libraries from:
1. **Zip file input** - Using a downloaded MSigDB zip file
2. **Direct download** - Downloading directly from MSigDB (requires login)

**NEW**: The scripts now also automatically update background gene lists with all unique genes from the C2 and C5 GMT files, validated against the Homo_sapiens.gene_info database.

**IMPORTANT**: These scripts are configured to process **only C2 and C5 gene sets** (excluding c2.all and c5.all) as these are the supported gene set categories for the enrichment tool.

## Files

- `update_msigdb.py` - Main update script with full functionality
- `update_msigdb_simple.py` - Simple wrapper script for the zip file in `data/recent_release/`
- `MSIGDB_UPDATE_README.md` - This documentation file

## Quick Start

If you have the MSigDB zip file in the `data/recent_release/` directory:

```bash
python3 update_msigdb_simple.py
```

This will:
1. Find the zip file at `data/recent_release/msigdb_v2025.1.Hs_files_to_download_locally.zip`
2. Update C2 and C5 gene set libraries (excluding c2.all and c5.all)
3. Download and process NCBI gene information files
4. Create a comprehensive background gene list
5. Update all necessary alias files

## Full Usage

### Using a Zip File

```bash
python update_msigdb.py --zip-file path/to/msigdb_v2025.1.Hs_files_to_download_locally.zip
```

### Direct Download (Requires Login)

```bash
python update_msigdb.py --download --email your.email@example.com --password your_password
```

### Options

- `--zip-file PATH` - Path to MSigDB zip file
- `--download` - Download directly from MSigDB
- `--email EMAIL` - Email for MSigDB login (required with --download)
- `--password PASSWORD` - Password for MSigDB login (required with --download)
- `--no-backup` - Skip creating backup of existing libraries
- `--skip-background-update` - Skip updating background gene lists

## What the Scripts Do

1. **Backup Creation**: Creates a timestamped backup of existing libraries in `data/libraries/backup/`
2. **File Extraction**: Extracts the MSigDB zip file to a temporary directory
3. **File Filtering**: Only processes C2 and C5 human symbols files (`*c2*.Hs.symbols.gmt` and `*c5*.Hs.symbols.gmt`, excluding c2.all and c5.all)
4. **Library Update**: Copies new/updated C2 and C5 .gmt files to `data/libraries/` (excluding c2.all and c5.all)
5. **Alias Update**: Updates `data/libraries/alias.json` with new library information
6. **NEW: Gene Collection**: Collects all unique genes from C2 and C5 GMT files (excluding c2.all and c5.all)
7. **NEW: Gene Validation**: Validates genes against Homo_sapiens.gene_info database
8. **NEW: Background Update**: Creates/updates background gene lists with validated genes
9. **Cleanup**: Removes temporary files

## File Structure

After running the update, your directories will contain:

### Libraries Directory
```
data/libraries/
├── alias.json                    # Library aliases and metadata
├── backup/                       # Backup of previous libraries
│   └── backup_YYYYMMDD_HHMMSS/
├── c2.cp.v2025.1.Hs.symbols.gmt  # C2: CP: Canonical pathways
├── c2.cp.biocarta.v2025.1.Hs.symbols.gmt # C2: CP: BioCarta
├── c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt # C2: CP: KEGG MEDICUS
├── c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt # C2: CP: KEGG Legacy
├── c2.cp.pid.v2025.1.Hs.symbols.gmt # C2: CP: PID
├── c2.cp.reactome.v2025.1.Hs.symbols.gmt # C2: CP: Reactome
├── c2.cp.wikipathways.v2025.1.Hs.symbols.gmt # C2: CP: WikiPathways
├── c2.cgp.v2025.1.Hs.symbols.gmt # C2: CGP: Chemical & genetic perturbations
├── c5.go.v2025.1.Hs.symbols.gmt  # C5: GO: Gene Ontology
├── c5.go.bp.v2025.1.Hs.symbols.gmt # C5: GO: Biological Process
├── c5.go.cc.v2025.1.Hs.symbols.gmt # C5: GO: Cellular Component
├── c5.go.mf.v2025.1.Hs.symbols.gmt # C5: GO: Molecular Function
├── c5.hpo.v2025.1.Hs.symbols.gmt # C5: HPO: Human Phenotype Ontology
└── ... (other C2 and C5 gene set files, excluding c2.all and c5.all)
```

### Backgrounds Directory
```
data/backgrounds/
├── alias.json                    # Background aliases and metadata
├── all_genes_from_c2_c5_msigdb_no_all.txt # All unique genes from C2 and C5 MSigDB (excluding c2.all and c5.all) (NEW)
├── all_genes_from_c2_c5_msigdb.txt # Previous C2/C5 background
├── all_genes_from_msigdb.txt     # Previous all-gene background
└── hgnc_symbols_2023-01-01.txt   # Original background file
```

## Background Gene List Update

### What's New

The scripts now automatically:

1. **Collect All Genes**: Extract all unique gene symbols from C2 and C5 GMT files (excluding c2.all and c5.all)
2. **Validate Genes**: Check each gene against the Homo_sapiens.gene_info database
3. **Create Background**: Generate a comprehensive background file with all validated genes
4. **Update Aliases**: Add the new background to the alias.json file

### Validation Process

- Genes are validated against the Homo_sapiens.gene_info database
- Only genes that exist in the database are included in the background
- Invalid genes are logged for review
- The process ensures compatibility with the gene enrichment analysis

### Background File Details

- **File**: `data/backgrounds/all_genes_from_c2_c5_msigdb_no_all.txt`
- **Format**: One gene symbol per line
- **Content**: All unique, validated genes from C2 and C5 MSigDB GMT files (excluding c2.all and c5.all)
- **Size**: Typically 20,000+ genes (varies by MSigDB version)

## Supported Gene Set Categories

The scripts are configured to process only the following MSigDB gene set categories (excluding c2.all and c5.all):

### C2: Curated Gene Sets
- **C2.cp**: Canonical pathways
- **C2.cp.biocarta**: BioCarta pathways
- **C2.cp.kegg_medicus**: KEGG MEDICUS pathways
- **C2.cp.kegg_legacy**: KEGG Legacy pathways
- **C2.cp.pid**: PID pathways
- **C2.cp.reactome**: Reactome pathways
- **C2.cp.wikipathways**: WikiPathways
- **C2.cgp**: Chemical & genetic perturbations

### C5: GO Gene Sets
- **C5.go**: Gene Ontology
- **C5.go.bp**: Biological Process
- **C5.go.cc**: Cellular Component
- **C5.go.mf**: Molecular Function
- **C5.hpo**: Human Phenotype Ontology

## Prerequisites

- Python 3.6+
- Required Python packages:
  - `requests` (for direct download functionality)

Install required packages:
```bash
pip install requests
```

## Troubleshooting

### Common Issues

1. **Zip file not found**
   - Ensure the zip file is in the correct location
   - Check file permissions

2. **Permission errors**
   - Ensure you have write permissions to the `data/libraries/` and `data/backgrounds/` directories

3. **Download fails**
   - Check your internet connection
   - Verify MSigDB login credentials
   - MSigDB may require registration for downloads

4. **No .gmt files found**
   - Ensure the zip file contains MSigDB gene set files
   - Check if the zip file is corrupted

5. **Gene validation issues**
   - Ensure Homo_sapiens.gene_info file exists in the project root
   - Check if the gene info file is corrupted or incomplete

6. **Background update fails**
   - Check if the backgrounds directory is writable
   - Verify that the gene converter module is accessible

### Logging

The scripts use Python's logging module. To see more detailed output, you can modify the logging level in the script:

```python
logging.basicConfig(level=logging.DEBUG)
```

## MSigDB Information

- **Website**: https://www.gsea-msigdb.org/
- **Current Version**: 2025.1.Hs
- **Download URL**: https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/msigdb_v2025.1.Hs_files_to_download_locally.zip

## Support

For issues with the update scripts, please check:
1. The logging output for error messages
2. File permissions and paths
3. MSigDB website status and login requirements
4. Homo_sapiens.gene_info file availability and integrity

For issues with the gene enrichment analysis application, refer to the main README.md file.
