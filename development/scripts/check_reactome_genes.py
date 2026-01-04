#!/usr/bin/env python3
"""
Script to check the unique number of genes in Reactome library.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.gene_set_library import GeneSetLibrary
from collections import Counter

def analyze_reactome_genes():
    """Analyze the Reactome library to show gene statistics."""
    
    print("🔬 REACTOME GENE ANALYSIS")
    print("=" * 50)
    
    # Load Reactome library
    try:
        reactome_lib = GeneSetLibrary(
            "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt",
            name="C2: CP: Reactome Pathways"
        )
    except Exception as e:
        print(f"Error loading Reactome library: {e}")
        return
    
    print(f"📚 Library: {reactome_lib.name}")
    print(f"📊 Total terms: {reactome_lib.num_terms}")
    print(f"🧬 Unique genes: {reactome_lib.size}")
    
    # Count how many terms each gene appears in
    gene_term_counts = Counter()
    for term in reactome_lib.library:
        for gene in term["genes"]:
            gene_term_counts[gene] += 1
    
    print(f"\n📈 Gene Distribution Statistics:")
    print(f"   Genes appearing in 1 term: {sum(1 for count in gene_term_counts.values() if count == 1)}")
    print(f"   Genes appearing in 2-5 terms: {sum(1 for count in gene_term_counts.values() if 2 <= count <= 5)}")
    print(f"   Genes appearing in 6-10 terms: {sum(1 for count in gene_term_counts.values() if 6 <= count <= 10)}")
    print(f"   Genes appearing in >10 terms: {sum(1 for count in gene_term_counts.values() if count > 10)}")
    
    # Show most common genes
    print(f"\n🏆 Top 10 most common genes (appearing in most terms):")
    for gene, count in gene_term_counts.most_common(10):
        print(f"   {gene}: {count} terms")
    
    # Show some term size statistics
    term_sizes = [term["size"] for term in reactome_lib.library]
    print(f"\n📏 Term Size Statistics:")
    print(f"   Smallest term: {min(term_sizes)} genes")
    print(f"   Largest term: {max(term_sizes)} genes")
    print(f"   Average term size: {sum(term_sizes) / len(term_sizes):.1f} genes")
    print(f"   Median term size: {sorted(term_sizes)[len(term_sizes)//2]} genes")
    
    # Verify the unique gene count
    all_genes = set()
    for term in reactome_lib.library:
        all_genes.update(term["genes"])
    
    print(f"\n✅ Verification:")
    print(f"   Manual count of unique genes: {len(all_genes)}")
    print(f"   Library's unique_genes size: {reactome_lib.size}")
    print(f"   Match: {'✅' if len(all_genes) == reactome_lib.size else '❌'}")

if __name__ == "__main__":
    analyze_reactome_genes()
