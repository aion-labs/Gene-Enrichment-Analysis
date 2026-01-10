# Streamlit Web Interface User Guide

This comprehensive guide explains how to use the Streamlit web interface for Iterative Gene Enrichment Analysis.

## Table of Contents

1. [Overview](#overview)
2. [Getting Started](#getting-started)
3. [Input Methods](#input-methods)
4. [Analysis Modes](#analysis-modes)
5. [Parameters and Settings](#parameters-and-settings)
6. [Understanding Results](#understanding-results)
7. [Export Options](#export-options)
8. [Advanced Features](#advanced-features)
9. [Troubleshooting](#troubleshooting)

## Overview

The Streamlit web interface provides an interactive, user-friendly way to perform gene enrichment analysis. It offers:

- **Interactive Web UI**: Access through any web browser
- **Real-time Visualization**: Charts, graphs, and tables update instantly
- **Multiple Input Methods**: Direct input, file selection, or file upload
- **Two Analysis Modes**: Regular and Iterative enrichment
- **Network Visualization**: Generate network graphs for AI systems
- **Statistical Benchmarking**: Compare results against random gene lists


**Recommendation:**
- **Use Streamlit** for interactive exploration, visualization, and one-off analyses
- **Use CLI** for batch processing, automation, and integration into computational pipelines

## Getting Started

### Launching the Application

#### Local Installation

```bash
# Activate virtual environment (if using one)
source .venv/bin/activate

# Run the Streamlit app
streamlit run code/streamlit_app.py
```

The application will automatically open in your default web browser at `http://localhost:8501`.

#### Docker Installation

```bash
# Build the Docker image
docker build -t cc_enrichment:latest .

# Run the container
docker run -p 8080:8080 cc_enrichment:latest
```

Access the application at `http://localhost:8080`.

### Interface Overview

The Streamlit interface consists of:

1. **Sidebar**: Mode selection and app information
2. **Main Panel**: 
   - **Analysis Tab**: Input area, settings, and parameters
   - **Advanced Settings Tab**: Additional configuration options
3. **Results Area**: Displays analysis results, charts, and download options

## Input Methods

### Method 1: Direct Input

1. **Select Input Format**: Choose between "Gene Symbols" or "Entrez IDs" from the dropdown
2. **Enter Genes**: Paste your gene list into the text area (one gene per line)
   - **Gene Symbols**: Official gene symbols (e.g., `TP53`, `BRCA1`, `EGFR`)
   - **Entrez IDs**: Numeric identifiers from NCBI Entrez (e.g., `7157`, `672`, `1956`)
3. **Name Your Gene Set**: Enter a descriptive name in the "Input/edit gene list name" field
4. **Click Submit**: Run the analysis

**Limitations:**
- Maximum 800 genes for optimal performance
- Minimum 20 genes required for analysis

### Method 2: File Selection

1. Use the dropdown menu: "Or select a file from the `data` folder"
2. Select a `.txt` file from the available gene lists
3. The file content will automatically populate the input area
4. The gene set name will be set to the filename (without extension)
5. Click Submit to run the analysis

**File Format:**
- Plain text files (`.txt`)
- One gene identifier per line
- Can contain gene symbols or Entrez IDs (must match selected format)

### Method 3: Load Example

1. Click the **"Load Example"** button
2. A pre-loaded example gene set will populate the input area
3. Modify if needed, then click Submit

### Gene Input Validation

The app automatically validates and converts your input:

- **Format Validation**: Checks gene identifiers against the selected format
- **Conversion**: 
  - Entrez IDs → Official gene symbols
  - Old symbols → Updated symbols (when possible)
- **Duplicate Removal**: Removes duplicate genes
- **Background Check**: Validates genes against the background list
- **Statistics Display**: Shows conversion statistics and unrecognized genes

**Validation Results Display:**
- ✅ **Validated genes**: Successfully converted and recognized
- ⚠️ **Unrecognized genes**: Genes that couldn't be matched (check spelling/format)
- 📊 **Statistics**: Total input, validated, converted, and unrecognized counts

## Analysis Modes

### Regular Mode

**Purpose**: Standard one-pass enrichment analysis that reports all significant results.

**How it works:**
This is a standard over-representation analysis (ORA)
1. Tests your gene set against all selected libraries
2. Calculates enrichment statistics for each term
3. Filters results based on p-value thresholds and overlap requirements
4. Displays ranked results with p-values and FDR adjustments

**Best for:**
- Standard enrichment analysis
- Getting a comprehensive view of all significant pathways
- When you want FDR-adjusted p-values

### Iterative Mode

**Purpose**: Iteratively removes genes from top hits to identify core enriched pathways.

**How it works:**
1. **Iteration 1**: Finds the top enriched term (lowest p-value)
2. **Removes genes**: Removes overlapping genes from your input gene set
3. **Iteration 2+**: Repeats enrichment with the reduced gene set
4. **Continues**: Until no terms pass the p-value threshold or max iterations reached
5. **Builds network**: Integrate enriched terms across library using a gene--term network.

**Best for:**
- Large language model assisted interpretation of enrichment analysis
- Visualizing enrichment results in a network

## Parameters and Settings

### Regular Mode Parameters

#### Raw P-value Threshold
- **Default**: `0.01`
- **Description**: Maximum raw p-value for terms to be included in results
- **Lower values**: More stringent filtering (fewer results)
- **Higher values**: Less stringent filtering (more results)

#### Adjusted P-value Threshold (FDR)
- **Default**: `0.05`
- **Description**: Maximum FDR-adjusted p-value for terms to be included
- **Note**: Only applicable in Regular mode

#### Minimum Overlap
- **Default**: `3`
- **Description**: Minimum number of genes that must overlap between your input and a term. The motivation here is to avoid inflating the list with many terms where only 1 or 2 genes overlap. Such terms might be significant statistically but typically have lower biological interpretation value.
- **Lower values**: Include terms with fewer overlapping genes
- **Higher values**: Require more overlap (more stringent)

#### Term Size Range
- **Default**: `10` to `600` genes
- **Description**: Filter gene sets by their size (number of genes). The motivation is to avoide removing genes with terms that are less informative for interpretation. Terms like "Receptors" can include 2000+ genes and their contribution to biological interpretation  would be low. 
- **Minimum**: Excludes very small gene sets
- **Maximum**: Excludes very large gene sets (reduces noise)

### Iterative Mode Parameters

#### P-value Threshold
- **Default**: `0.01`
- **Description**: Maximum p-value for terms to be selected at each iteration. Higher values are not supported for statistical benchmarking.
- **Note**: No FDR adjustment in iterative mode

#### Max Iterations
- **Default**: `10`
- **Description**: Maximum number of iterations to perform. The main motivation is saving compute time. Statistical analysis is supported up to 30 iteration.
- **`0`**: No limit (continues until no terms pass threshold)
- **Higher values**: Allows more iterations (longer analysis)


### Advanced Settings

#### Number of Results to Display
- **Default**: `10`
- **Description**: Controls how many top-ranked results are displayed in the interactive table and bar chart for each library. 
- **What it affects**:
  - ✅ **Table view**: Limits the number of rows shown in the results table
  - ✅ **Bar chart**: Limits the number of bars displayed in the visualization
  - ❌ **Download files**: Does NOT affect downloadable TSV/JSON files (these contain ALL results that pass the filters)
- **Use cases**:
  - **Small value (5-10)**: Focus on the most significant results, cleaner visualization
  - **Large value (50-100)**: See more results at once, useful for comprehensive exploration
- **Note**: Only applies to Regular mode. Iterative mode shows all iterations by default.

#### P-value Calculation Method
- **Options**:
  - `Fisher's Exact Test` (default, recommended)
  - `Hypergeometric Test`
  - `Chi-squared Test`
- **Description**: Statistical method used to calculate enrichment p-values
- **Recommendation**: Use Fisher's Exact Test unless you have a specific reason

#### Background Gene List
- **Default**: `all_genes.txt`
- **Description**: The reference gene set used for enrichment calculations
- **Options**: Select from available backgrounds in the dropdown
- **Custom Upload**: Upload your own background file (see Advanced Features)

#### Gene Set Libraries
- **Description**: Select which gene set libraries to test against
- **Selection**: Checkboxes for each available library
- **Default Libraries**:
  - H: Hallmark Gene Sets
  - C2: Reactome Pathways
  - C5: Gene Ontology: Biological Process
- **Available Libraries**: MSigDB collections (H, C2, C5, etc.)

## Understanding Results

### Regular Mode Results

#### Results Table
Each library displays a table with:

| Column | Description |
|--------|-------------|
| **Rank** | Ranking by p-value (1 = most significant) |
| **Term** | Gene set term name |
| **Description** | Detailed description of the term |
| **Overlap size** | Format: `X/Y` (overlapping genes / term size) |
| **Genes** | List of overlapping genes (click to expand) |
| **p-value** | Raw p-value from statistical test |
| **-log(p-value)** | Negative log10 of p-value (higher = more significant) |
| **FDR** | False Discovery Rate adjusted p-value |

#### Bar Chart
- **X-axis**: `-log10(p-value)` (higher = more significant)
- **Y-axis**: Term names
- **Interactive**: Hover for details, zoom, pan
- **Color**: Gradient based on significance

#### Library-Specific Information
Each library result shows:
- **Gene Count**: `X/Y genes` (library-specific input / total input)
- **Background Count**: `X/Y background` (library-specific background / total background)

### Iterative Mode Results

#### Iteration Table
Each library displays results organized by iteration:

| Column | Description |
|--------|-------------|
| **Iteration** | Iteration number (1, 2, 3, ...) |
| **Term** | Top enriched term for this iteration |
| **Description** | Term description |
| **iteration overlapping genes** | Format: `X/Y` (overlap / term size) |
| **iteration p-value** | P-value for this iteration |
| **iteration -log(p-value)** | Negative log10 of iteration p-value |
| **Genes removed for next iteration** | Genes removed from input for next iteration |
| **Full list overlapping genes** | Overlap with original (full) gene list |
| **Full list p-value** | P-value against full gene list |
| **Regular FDR** | FDR from iteration 1 (ORA comparison) |

#### Network Visualization
- **Interactive Graph**: Shows relationships between iterations
- **Nodes**: Represent enriched terms
- **Edges**: Show gene removal relationships
- **Layout**: Hierarchical structure showing pathway relationships

#### ORA vs iGEA Comparison

- **ORA**: Over-Representation Analysis (Iteration 1 only)
- **iGEA**: Iterative Gene Enrichment Analysis (All iterations)
- **Statistics**: Comparison metrics showing how iterative analysis differs from standard ORA


## Export Options

### Regular Mode Exports

#### Combined Results
- **TSV File**: All libraries combined in one tab-separated file
- **JSON File**: All libraries combined in JSON format

#### Individual Library Files
- **Archive**: Single `.tar.gz` file containing:
  - Individual TSV files for each library
  - Individual JSON files for each library

### Iterative Mode Exports

#### ORA Results (Iteration 1)
- **Combined TSV**: All libraries' iteration 1 results combined
- **Individual Files Archive**: All libraries' ORA files in one archive

#### iGEA Results (All Iterations)
- **Combined TSV**: All libraries' all-iteration results combined
- **Individual Files Archive**: All libraries' iteration files in one archive

#### Validated Gene Symbols
- **TXT File**: Your input genes after validation and conversion
- **Format**: One gene symbol per line

### File Formats

#### TSV (Tab-Separated Values)
- **Use**: Spreadsheet applications (Excel, Google Sheets)
- **Format**: Tab-separated columns
- **Encoding**: UTF-8

#### JSON (JavaScript Object Notation)
- **Use**: Programmatic access, data exchange
- **Format**: Structured JSON with metadata
- **Includes**: Full analysis parameters and results

#### DOT (Graph Description Language)
- **Use**: Network visualization, graph analysis
- **Format**: Graphviz DOT format
- **Available**: Only in Iterative mode (network generation)

## Advanced Features

### Custom Background Gene List

1. Go to **Advanced Settings** tab
2. Select background format: "Gene Symbols" or "Entrez IDs"
3. Click **"Upload your background gene list"**
4. Select a `.txt` file (one gene per line)
5. Click **"Apply Settings"**
6. Your custom background will appear in the dropdown menu

**Use Cases:**
- Tissue-specific backgrounds
- Experiment-specific gene sets
- Custom reference genomes

### Network Generation (Iterative Mode Only)

#### Purpose
Generate network graphs for conversational AI systems (e.g., ChatGPT, Claude).

#### Steps

1. **Run Iterative Analysis**: Complete an iterative enrichment analysis
2. **Select Libraries**: 
   - Check "Use results in network" for each library you want to include
   - Or use the "Generate network for conversational AI systems" section
3. **Generate Network**: Click "Generate Network" button
4. **View Network**: Interactive Plotly graph displays
5. **AI Prompt**: Copy the generated prompt for use with AI systems

#### Network Features
- **Interactive Visualization**: Zoom, pan, hover for details
- **Hierarchical Layout**: Shows pathway relationships
- **Gene Count Display**: Shows number of genes in your input
- **Library List**: Lists all included libraries
- **DOT Format**: Available for download

### Statistical Benchmarking (Iterative Mode Only)

#### Purpose
Compare your network connectivity against random gene lists of similar size.

#### How It Works

1. **Automatic Computation**: Runs automatically after iterative analysis
2. **Null Distribution**: Uses pre-computed permutation data
3. **Comparison**: Compares your network metrics against random gene lists
4. **Statistics**: Provides Z-scores and percentiles

#### Metrics Compared

| Metric | Description |
|--------|-------------|
| **Genes** | Number of unique genes in network |
| **Terms** | Number of enriched terms |
| **Edges** | Number of connections in network |
| **Density** | Network connectivity density |

#### Interpreting Results

**Z-score and Percentile:**
- **Z-score > 2.0 and Percentile > 95%**: ✓✓ SIGNIFICANTLY BETTER (top 5%)
- **Z-score > 1.0 and Percentile > 84%**: ✓ BETTER (top 16%)
- **Z-score ~ 0.0 and Percentile ~ 50%**: ~ SIMILAR to random
- **Z-score < -1.0 and Percentile < 16%**: ✗ WORSE than random

**What This Means:**
- Higher values indicate better network connectivity
- Your gene list is compared against thousands of random gene lists
- Significant results suggest biological coherence

#### Requirements

- **Permutation Data**: Requires pre-computed permutation statistics
- **Gene List Size**: Must match available permutation data sizes
- **Parameters**: Best results with default parameters (min overlap: 3, term size: 10-600)

#### Download Statistical Report

- **Format**: Plain text file
- **Content**: Comprehensive statistical analysis report
- **Includes**: All benchmark metrics, interpretations, and library information

### MSigDB Version Information

The app displays:
- **MSigDB Version**: Current version of gene set libraries (e.g., v2025.1)
- **Last Update Date**: When libraries were last updated

This information appears at the top of the Analysis tab.

## Troubleshooting

### Common Issues

#### "No enrichment results found"

**Possible Causes:**
1. **P-value threshold too strict**: Try increasing the p-value threshold
2. **Minimum overlap too high**: Reduce the minimum overlap requirement
3. **Term size filters too restrictive**: Adjust min/max term size range
4. **No significant enrichment**: Your gene list may not be enriched in selected libraries

**Solutions:**
- Increase p-value threshold to `0.05` or higher
- Reduce minimum overlap to `1` or `2`
- Expand term size range (e.g., `5` to `1000`)
- Try different gene set libraries
- Check that your genes are valid and in the background

#### "Gene validation failed"

**Possible Causes:**
1. **Invalid gene identifiers**: Genes not recognized in database
2. **Wrong format selected**: Using Entrez IDs but selected "Gene Symbols" (or vice versa)
3. **Too few genes**: Less than 20 validated genes
4. **Too many genes**: More than 800 genes

**Solutions:**
- Check gene spelling and format
- Verify you selected the correct input format
- Add more genes (minimum 20 required)
- Reduce gene list size (maximum 500)
- Check unrecognized genes list and correct them

#### "Background validation failed"

**Possible Causes:**
1. **Background file missing**: Background file not found
2. **Background file corrupted**: File cannot be read
3. **Format mismatch**: Background format doesn't match selection

**Solutions:**
- Refresh the page
- Try selecting a different background
- Check that background files exist in `data/backgrounds/`
- Re-upload custom background if using one

#### "Library validation failed"

**Possible Causes:**
1. **No libraries selected**: No gene set libraries checked
2. **Library file missing**: Selected library file not found
3. **Library file corrupted**: Library file cannot be read

**Solutions:**
- Select at least one library from the checkboxes
- Refresh the page
- Check that library files exist in `data/libraries/`

#### App is slow or unresponsive

**Possible Causes:**
1. **Large gene list**: Too many genes (>800)
2. **Many libraries selected**: Processing many libraries simultaneously
3. **Large term sizes**: Very large gene sets in libraries
4. **Browser issues**: Browser memory or performance

**Solutions:**
- Reduce gene list size
- Select fewer libraries
- Adjust max term size filter
- Refresh the browser
- Close other browser tabs
- Use CLI for large-scale analyses

#### Network generation fails

**Possible Causes:**
1. **No libraries selected**: No libraries chosen for network
2. **No iterative results**: Must run iterative analysis first
3. **No enrichment results**: Libraries have no significant terms

**Solutions:**
- Run iterative analysis first
- Select at least one library with results
- Check that libraries have enrichment results
- Adjust parameters to get more results

#### Benchmarking not available

**Possible Causes:**
1. **No permutation data**: Permutation statistics files missing
2. **Gene list size mismatch**: Your gene list size not in permutation data
3. **Parameter mismatch**: Using non-default parameters
4. **No clusters**: Network has no clusters to benchmark

**Solutions:**
- Check that permutation data exists in `permutations/` directory
- Use default parameters for best compatibility
- Ensure gene list size matches available permutation data
- Verify network has at least one cluster

### Performance Tips

1. **Start Small**: Begin with a small gene list and few libraries
2. **Use Filters**: Apply term size filters to reduce computation
3. **Close Unused Tabs**: Free up browser memory
4. **Regular Mode First**: Use Regular mode to quickly test parameters
5. **Iterative Mode**: Use Iterative mode only when needed (more computationally intensive)

### Getting Help

If you encounter issues not covered here:

1. **Check Logs**: Look at terminal/console output for error messages
2. **Validate Input**: Ensure gene list format is correct
3. **Try Defaults**: Reset to default parameters
4. **Check Documentation**: Review CLI documentation for parameter details
5. **File Issues**: Report bugs with error messages and input examples

## Best Practices

### Input Preparation

1. **Use Official Gene Symbols**: Prefer official HGNC gene symbols
2. **Check Format**: Verify Entrez IDs are numeric
3. **Remove Duplicates**: Clean your gene list before input
4. **Reasonable Size**: Keep gene lists between 20-800 genes
5. **Meaningful Names**: Use descriptive gene set names

### Analysis Strategy

1. **Start with Regular Mode**: Get overview of all significant pathways
2. **Use Default Libraries**: Start with Hallmark, Reactome, GO BP
3. **Adjust Parameters Gradually**: Make small parameter changes
4. **Iterative for Networks**: Use Iterative mode when you need networks
5. **Save Results**: Download results for record-keeping

### Result Interpretation

1. **Check P-values**: Focus on highly significant results (p < 0.01)
2. **Consider FDR**: In Regular mode, FDR < 0.05 is significant
3. **Review Overlap**: Check that overlap sizes are meaningful
4. **Compare Libraries**: Results may vary between libraries
5. **Biological Context**: Interpret results in context of your research

### Export and Documentation

1. **Save TSV Files**: For spreadsheet analysis
2. **Keep JSON Files**: For programmatic access
3. **Document Parameters**: Note parameter settings used
4. **Record Gene Lists**: Save validated gene symbols
5. **Version Control**: Note MSigDB version used

---

## Quick Reference

### Keyboard Shortcuts
- **Ctrl/Cmd + R**: Refresh page (resets app)
- **Ctrl/Cmd + F**: Find in page

### Default Settings Summary

**Regular Mode:**
- Raw p-value: `0.01`
- Adjusted p-value (FDR): `0.05`
- Min overlap: `3`
- Term size: `10-600`
- Results to display: `10`

**Iterative Mode:**
- P-value: `0.01`
- Max iterations: `10`
- Min overlap: `3`
- Term size: `10-600`

**Advanced:**
- P-value method: `Fisher's Exact Test`
- Background: `all_genes.txt`
- Input format: `Gene Symbols`

### File Locations

- **Results**: `results/run_YYYYMMDD_HHMMSS/`
- **Gene Lists**: `data/gene_lists/`
- **Backgrounds**: `data/backgrounds/`
- **Libraries**: `data/libraries/`

---

**Last Updated**: Based on app version with MSigDB v2025.1 support

For CLI documentation, see [CLI_README.md](CLI_README.md)
