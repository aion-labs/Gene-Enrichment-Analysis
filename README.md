# Iterative Gene Enrichment Analysis

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://iterative-gene-enrichment.streamlit.app/)

> **Try it now:** A slim hosted version of this tool is available at [iterative-gene-enrichment.streamlit.app](https://iterative-gene-enrichment.streamlit.app/) — no installation required.

## Prerequisites
- Have a list of genes for enrichment analysis in official human gene symbols or entrez IDs.
- Unless you use all genes as a background, have a background file with all measured genes ready in the human gene symbols or entrez IDs.
- You may add your own gene set collections in GMT format with `python code/libraries.py add` (see "Managing gene-set collections").
- NCBI gene data used to validate symbols ship with the repository (`data/recent_release/Homo_sapiens.gene_info.gz`,
  `gene_history_9606.gz`). iGEA stops with an error if they are missing, so results never depend on the working directory.
- One universe: input genes outside the background are removed (and listed), and gene sets are restricted to the
  background before the size filter (10-600 genes), so the input, the gene sets and the random lists share the background.

## Run the Application

### Web Interface (Streamlit)
The application provides a user-friendly web interface that runs in your browser:
- **To run locally**: `streamlit run code/streamlit_app.py`
- **To run in Docker**: `docker build -t cc_enrichment:latest .` and `docker run -p 8080:8080 cc_enrichment:latest`

The Streamlit interface offers:
- Interactive web-based UI accessible through your browser
- Real-time parameter adjustment
- Visual results display with charts and graphs
- Easy file upload and download capabilities

### Command Line Interface (CLI)
For automated analysis, batch processing, or integration into workflows, use the command-line interface:

```bash
# Basic usage
python code/cli.py --genelist your_genes.txt --libraries "H: Hallmark Gene Sets"

# Run with all active libraries
python code/cli.py --genelist your_genes.txt

# Iterative enrichment mode
python code/cli.py --genelist your_genes.txt --mode iterative

# See all options
python code/cli.py --help
```

**CLI Features:**
- Full feature parity with the web interface
- Support for both regular and iterative enrichment modes
- Batch processing capabilities
- Automated output generation (TSV, JSON, DOT files)
- Gene format support (symbols or Entrez IDs)
- Configurable statistical parameters

For detailed CLI documentation, see [CLI_README.md](documentation/CLI_README.md).

### Choosing Between Web Interface and CLI

| Feature | Web Interface (Streamlit) | Command Line Interface (CLI) |
|---------|---------------------------|------------------------------|
| **Ease of Use** | ✅ User-friendly, no command line knowledge required | ⚠️ Requires basic command line skills |
| **Visualization** | ✅ Interactive charts and graphs | ❌ Text-based output only |
| **Real-time Adjustments** | ✅ Instant parameter changes | ❌ Requires re-running commands |
| **Batch Processing** | ❌ Manual operation only | ✅ Automated processing |
| **Integration** | ❌ Limited to web browser | ✅ Easy integration into workflows |
| **Resource Usage** | ⚠️ Higher memory usage | ✅ Lower resource footprint |
| **Remote Access** | ✅ Accessible from any browser | ⚠️ Requires SSH/remote access |

**Recommendation:**
- **Use Streamlit** for interactive exploration, visualization, and one-off analyses
- **Use CLI** for batch processing, automation, and integration into computational pipelines

## Statistics (per-run random gene lists)

iGEA evaluates each analysis against random gene lists generated **for that run**: each random list has
exactly as many genes as your validated input list, is drawn without replacement from **your**
background, and is analysed with the same collections, test, p-value threshold, minimum overlap,
gene-set size limits and iteration cap. Two layers of statistics are applied, identically to your list
and to every random list:

1. **Term level — signal-adjusted empirical q-value.** For a selected term with p-value p, q is the
   expected number of chance terms with p' ≤ p (from the random lists, per collection) divided by the
   number of your terms with p' ≤ p, made monotone. Genes assigned to accepted terms are then treated as
   signal and the expected chance terms are re-estimated from random lists of the remaining number of
   genes, until stable. **The network is built only from terms with q < 0.10 (default).** Use a q-value
   threshold of 1 to keep all selected terms (q-values are still reported).
2. **Network level.** The largest connected cluster of the filtered network is compared with the
   largest clusters of the filtered random networks: z-score and empirical p = (r + 1) / (n + 1), where r
   is the number of random lists whose cluster is at least as large and n the number of random lists
   (default 10,000; smallest attainable p ≈ 0.0001).

Random lists are reproducible from the seed (default 2026). A fast vectorised engine
(`code/fast_igea.py`, identical results to the reference implementation) makes this practical: for a
272-gene list, 10,000 random lists plus the q-value adjustment take about 1.5 min on 4 cores.

CLI: `python code/cli.py -g genes.txt -m iterative [--null-lists 10000] [--q-threshold 0.10] [--null-seed 2026] [--null-jobs 0] [--null-cache DIR]`.
`--null-cache DIR` keeps the generated random lists, so a run interrupted by a time limit resumes when the same
command is rerun (identical results). Iterations stop when no gene set reaches p < 0.01 or after 30 iterations
per collection (`--max-iterations`, 0 = no limit). Lists above 800 genes run with a warning (longer run time).
Statistics run by default in iterative mode; `--q-threshold 1` keeps all selected terms; `--no-benchmark` skips
statistics (network from all terms, unfiltered); `--null-source stored` reproduces the previous behaviour
(precomputed permutations, default background only).

## Managing gene-set collections

Collections live in `data/libraries/` and are listed in `data/libraries/alias.json` (name, file, active,
colour). Active collections are used by default by the CLI and the app (the 10 used in the paper are active).

```
python code/libraries.py list                                        # * = active
python code/libraries.py activate "C5: Human Phenotype Ontology"
python code/libraries.py deactivate "C2: BioCarta"
python code/libraries.py add my_sets.gmt --name "Custom: My sets"    # --inactive, --color "#1F77B4"
python code/libraries.py remove "Custom: My sets"                    # file moved to data/libraries/previous/removed/
```

Updating to a new MSigDB release (download the release zip from the MSigDB website, free login) and,
first, the NCBI gene data it should be harmonised with (`Homo_sapiens.gene_info.gz` and `gene_history.gz`
from NCBI's `gene/DATA` area):

```
python code/libraries.py update-gene-data --gene-info Homo_sapiens.gene_info.gz --gene-history gene_history.gz
python code/libraries.py update-msigdb --zip msigdb_v2026.1.Hs_files_to_download_locally.zip --dry-run
python code/libraries.py update-msigdb --zip msigdb_v2026.1.Hs_files_to_download_locally.zip
```

Every MSigDB collection listed in `alias.json` is replaced by its new version (`--only c2.cp.reactome ...` or
`--active-only` to limit), keeping its name, colour and active flag; gene symbols are mapped to current NCBI
symbols (unmapped symbols are kept); previous files are moved to `data/libraries/previous/vYYYY.N/`. The
background (`data/backgrounds/all_genes.txt`) is not changed. Results depend on the MSigDB release and NCBI
gene data, so record both; the paper uses MSigDB 2025.1 with the bundled gene data (tag `csbj-v7-revision`).

## Communities, hubs and LLM interpretation

Each iterative run also writes:
- `network_modules_*.tsv`: communities of the network (Leiden, modularity, on terms linked by shared genes), gene
  degree (= number of collections placing the gene in the network), participation coefficient and bridge genes.
- `combined_network.dot`: the network for drawing (colours by collection).
- `network_llm.dot`: the same network for reading (no styling): one cluster per community, term and gene
  attributes (library, iteration, p, q, community, degree, bridge) and the statistics as graph attributes.
- `llm_prompt.txt`: the prompt for a language model: statistics (with an instruction for networks that are not
  distinguishable from random lists), community table, the semantic DOT, and instructions requiring every
  statement to cite network terms and genes, with outside knowledge marked [outside the network].

To send the prompt to an OpenAI model (standard library only; the key stays in your environment):
```
export OPENAI_API_KEY="sk-..."
python3 code/llm_client.py <run folder> --model <model-id> [--repeats 3]   # or: cli.py ... --llm-model <model-id>
```
Responses, full request/response JSON and a manifest (model, response id, time, tokens) are saved in `<run folder>/llm/`.

## Overview of Application Components
### Key Features

#### Gene Set Submission
Users have the flexibility to input their gene set in multiple ways:
- **Direct Input**: Paste a newline-separated list of gene identifiers and assign a name to the gene set.
- **File Selection**: Choose an existing gene list from the local `data` folder.
- **Multiple Formats**: Support for gene symbols (e.g., TP53, BRCA1) and Entrez IDs (e.g., 7157, 672).
- **Format Selection**: Users must specify their input format before entering gene lists.

#### Gene Set Validation
The app performs crucial validation steps to ensure the accuracy of the analysis:
- **Format Validation**: Validates gene identifiers against the selected input format (symbols or Entrez IDs).
- **Conversion**: Converts Entrez IDs to official gene symbols when using Entrez ID format. In addition, converts old symbols to updated symbols when possible. 
- **Duplicate Check**: Checks for any duplicate gene names within the input gene set.
- **Background Check**: Validates the input gene set against a predefined background gene list, which helps in identifying any genes that might not be present in the background list.

#### Results Display and Interaction
- Interactive Bar Graph: An interactive bar graph provides a visual representation of the results for each gene set library.

#### Data Export Options
- Individual Library Export: Users can save the results for each gene set library in TSV and JSON formats.
- Consolidated Export: There is an option to save results for all libraries in a single TSV file.

#### Results Management
All results, along with the complementary metadata, are stored in the `results` folder and are assigned unique filenames.

### Advanced Features

#### Customization of Results Display
Users can customize the number of results to display in both the chart and the bar graph in the Streamlit app.

#### P-value Calculation Methods
The app offers various statistical methods to calculate p-values:
- Fisher's Exact Test.
- Hypergeometric Test.
- Chi-Squared Test.

**Important Note on P-values**: The p-value threshold used for filtering results are both **raw p-value** and FDR adjusted p-values. In Iterative mode, adjusted p-values are not relevant and thus not calculated. 

#### Custom Reference Data Upload
Background Gene List Upload: Users can upload their background gene lists as .txt files containing newline-separated gene names or Entrez IDs.
Gene Set Libraries Upload: The app allows for the upload of gene set libraries in the .gmt file format, enabling the use of custom gene sets in the analysis.
