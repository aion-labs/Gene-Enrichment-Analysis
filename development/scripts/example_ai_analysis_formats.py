#!/usr/bin/env python3
"""
Example script to demonstrate the AI analysis formats for iterative enrichment networks.
This shows how the different formats look with sample data.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.ui.rendering import generate_ai_analysis_prompt, generate_structured_network_analysis, generate_json_network_analysis

def create_sample_dot():
    """Create a sample DOT network for demonstration."""
    return '''graph iterative_enrichment {
  graph [layout=neato];
  node [shape=ellipse];
  "gene_TP53" [label="TP53", type="gene"];
  "gene_BRCA1" [label="BRCA1", type="gene"];
  "gene_CDKN1A" [label="CDKN1A", type="gene"];
  "gene_BAX" [label="BAX", type="gene"];
  "gene_CASP3" [label="CASP3", type="gene"];
  "gene_ATM" [label="ATM", type="gene"];
  "gene_CHEK2" [label="CHEK2", type="gene"];
  "gene_RAD51" [label="RAD51", type="gene"];
  "term_1_GOBP_APOPTOSIS" [label="GOBP: APOPTOSIS", style=filled, fontcolor="white", fillcolor="#FF6B6B", type="term"];
  "term_2_REACTOME_DNA_REPAIR" [label="REACTOME: DNA REPAIR", style=filled, fontcolor="white", fillcolor="#4ECDC4", type="term"];
  "term_3_KEGG_CELL_CYCLE" [label="KEGG: CELL CYCLE", style=filled, fontcolor="white", fillcolor="#45B7D1", type="term"];
  "term_4_HALLMARK_STRESS_RESPONSE" [label="HALLMARK: STRESS RESPONSE", style=filled, fontcolor="white", fillcolor="#96CEB4", type="term"];
  "gene_TP53" -- "term_1_GOBP_APOPTOSIS";
  "gene_TP53" -- "term_2_REACTOME_DNA_REPAIR";
  "gene_TP53" -- "term_3_KEGG_CELL_CYCLE";
  "gene_TP53" -- "term_4_HALLMARK_STRESS_RESPONSE";
  "gene_BRCA1" -- "term_2_REACTOME_DNA_REPAIR";
  "gene_BRCA1" -- "term_3_KEGG_CELL_CYCLE";
  "gene_CDKN1A" -- "term_1_GOBP_APOPTOSIS";
  "gene_CDKN1A" -- "term_3_KEGG_CELL_CYCLE";
  "gene_BAX" -- "term_1_GOBP_APOPTOSIS";
  "gene_CASP3" -- "term_1_GOBP_APOPTOSIS";
  "gene_ATM" -- "term_2_REACTOME_DNA_REPAIR";
  "gene_ATM" -- "term_4_HALLMARK_STRESS_RESPONSE";
  "gene_CHEK2" -- "term_2_REACTOME_DNA_REPAIR";
  "gene_CHEK2" -- "term_4_HALLMARK_STRESS_RESPONSE";
  "gene_RAD51" -- "term_2_REACTOME_DNA_REPAIR";
}'''

def main():
    """Generate example AI analysis formats."""
    print("🎯 Generating Example AI Analysis Formats for Iterative Enrichment Networks")
    print("=" * 80)
    
    # Create sample DOT data
    sample_dot = create_sample_dot()
    
    print("\n📊 SAMPLE NETWORK DATA:")
    print("This example shows a network with 8 genes and 4 biological terms from 4 iterations")
    print("Genes: TP53, BRCA1, CDKN1A, BAX, CASP3, ATM, CHEK2, RAD51")
    print("Terms: GOBP: APOPTOSIS (iter 1), REACTOME: DNA REPAIR (iter 2), KEGG: CELL CYCLE (iter 3), HALLMARK: STRESS RESPONSE (iter 4)")
    print("\n" + "=" * 80)
    
    # Generate Enhanced Prompt Format
    print("\n🔬 FORMAT 1: ENHANCED PROMPT FORMAT")
    print("-" * 50)
    enhanced_prompt = generate_ai_analysis_prompt(sample_dot)
    print(enhanced_prompt[:1000] + "...\n[truncated for display]")
    
    # Generate Structured Table Format
    print("\n📋 FORMAT 2: STRUCTURED TABLE FORMAT")
    print("-" * 50)
    structured_format = generate_structured_network_analysis(sample_dot)
    print(structured_format[:1000] + "...\n[truncated for display]")
    
    # Generate JSON Format
    print("\n🔧 FORMAT 3: JSON FORMAT")
    print("-" * 50)
    json_format = generate_json_network_analysis(sample_dot)
    print(json_format[:1000] + "...\n[truncated for display]")
    
    print("\n" + "=" * 80)
    print("📁 To see the complete formats, the files would be saved as:")
    print("   • ai_analysis_prompt_enhanced.txt")
    print("   • ai_analysis_structured.txt") 
    print("   • ai_analysis_json.txt")
    print("\n💡 Each format provides the same network data but optimized for different AI analysis approaches!")

if __name__ == "__main__":
    main()
