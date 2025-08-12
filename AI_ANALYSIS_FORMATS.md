# AI Analysis Formats for Iterative Enrichment Networks

This document describes the different formats available for providing iterative enrichment network data to AI tools like ChatGPT for biological interpretation.

## Overview

The iterative enrichment analysis generates network data that can be exported in multiple formats optimized for different AI analysis approaches. Each format provides the same underlying data but structures it differently for optimal AI interpretation.

## Available Formats

### 1. Enhanced Prompt Format (`ai_analysis_prompt_enhanced.txt`)

**Best for**: Comprehensive biological analysis with detailed instructions

**Features**:
- Pre-processed network statistics and centrality analysis
- Structured biological interpretation guidelines
- Clear analysis request with specific questions
- Includes original DOT format for reference

**Content Structure**:
- Experimental context and methodology explanation
- Network statistics (genes, terms, connections, iterations)
- Gene centrality analysis (top 10 most connected genes)
- Iteration sequence analysis with significance levels
- Detailed analysis instructions with 6 specific sections
- Original DOT network data

**Use Case**: When you want ChatGPT to provide a comprehensive biological interpretation with specific insights about hub genes, pathway relationships, and clinical relevance.

### 2. Structured Table Format (`ai_analysis_structured.txt`)

**Best for**: Data-focused analysis with clear tabular organization

**Features**:
- Tabular data presentation for easy parsing
- Gene centrality rankings with connection counts
- Iteration sequence tables with significance levels
- Connection matrices showing gene-term relationships
- Biological terms summary with gene lists

**Content Structure**:
- Network summary statistics
- Gene centrality table (ranked by connections)
- Iteration sequence table (by significance)
- Connection matrix (gene → biological terms)
- Biological terms summary (with connected genes)
- Analysis instructions

**Use Case**: When you prefer a more data-driven approach where the AI can easily parse specific relationships and rankings.

### 3. JSON Format (`ai_analysis_json.txt`)

**Best for**: Machine-readable structured data analysis

**Features**:
- Fully structured JSON data format
- Hierarchical organization of network components
- Easy programmatic parsing
- Complete network topology information
- Hub gene analysis with biological roles

**Content Structure**:
- Network summary with key metrics
- Gene centrality array with rankings and connected terms
- Iteration sequence with significance levels
- Biological terms with gene counts and connections
- Hub genes analysis with centrality scores
- Analysis instructions

**Use Case**: When you want the most structured format for AI analysis or when building automated analysis pipelines.

## Format Comparison

| Feature | Enhanced Prompt | Structured Table | JSON Format |
|---------|----------------|------------------|-------------|
| **Ease of Reading** | High | Medium | Low (requires JSON parsing) |
| **Data Structure** | Narrative + Data | Tabular | Hierarchical JSON |
| **AI Parsing** | Good | Excellent | Best |
| **Biological Context** | High | Medium | Medium |
| **Statistical Analysis** | Included | Included | Included |
| **Network Topology** | Summary | Detailed | Complete |
| **File Size** | Medium | Small | Medium |

## Choosing the Right Format

### Use Enhanced Prompt Format when:
- You want comprehensive biological interpretation
- You prefer narrative-style analysis
- You need detailed instructions for the AI
- You want to include biological context and methodology

### Use Structured Table Format when:
- You prefer data-driven analysis
- You want to focus on specific relationships
- You need clear rankings and metrics
- You want easily scannable information

### Use JSON Format when:
- You're building automated analysis tools
- You need programmatic access to network data
- You want the most structured representation
- You plan to integrate with other computational tools

## Example Usage with ChatGPT

### Enhanced Prompt Format:
```
Copy the content of ai_analysis_prompt_enhanced.txt and paste it into ChatGPT with:
"Please analyze this iterative enrichment network and provide the requested biological interpretation."
```

### Structured Table Format:
```
Copy the content of ai_analysis_structured.txt and paste it into ChatGPT with:
"Please analyze this structured network data and provide insights about the biological processes and gene relationships."
```

### JSON Format:
```
Copy the content of ai_analysis_json.txt and paste it into ChatGPT with:
"Please analyze this JSON network data and provide a comprehensive biological interpretation focusing on hub genes and pathway relationships."
```

## Key Biological Insights to Look For

Regardless of format, the AI should identify:

1. **Hub Genes**: Genes with high connectivity that may be key regulators
2. **Biological Themes**: Common pathways or processes across iterations
3. **Iteration Patterns**: How biological significance changes across iterations
4. **Network Architecture**: Clusters, modules, and connectivity patterns
5. **Clinical Relevance**: Potential therapeutic targets or biomarkers
6. **Experimental Context**: Hypothesized experimental conditions

## Technical Notes

- All formats are generated from the same DOT network representation
- Statistical calculations (centrality, connectivity) are consistent across formats
- The original DOT format is preserved for visualization purposes
- All formats include the same core network data with different presentations

## Future Enhancements

Potential improvements to consider:
- Integration with external databases for gene function annotation
- Pathway enrichment analysis within the network
- Comparative analysis between different experimental conditions
- Machine learning-based network analysis recommendations
