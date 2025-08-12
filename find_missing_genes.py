#!/usr/bin/env python3
"""
Script to find genes that are in Reactome but not in the background of all genes.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.gene_set_library import GeneSetLibrary
from code.background_gene_set import BackgroundGeneSet

def find_missing_genes():
    """Find genes that are in Reactome but not in the background."""
    
    print("🔍 FINDING MISSING GENES")
    print("=" * 50)
    
    # Load Reactome library
    print("📚 Loading Reactome library...")
    try:
        reactome_lib = GeneSetLibrary(
            "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt",
            name="C2: CP: Reactome Pathways"
        )
        print(f"Reactome unique genes: {reactome_lib.size}")
    except Exception as e:
        print(f"Error loading Reactome library: {e}")
        return
    
    # Load background
    print("\n📊 Loading background genes...")
    try:
        background = BackgroundGeneSet("data/backgrounds/all_genes.txt", name="All Genes")
        print(f"Background genes: {background.size}")
    except Exception as e:
        print(f"Error loading background: {e}")
        return
    
    # Find genes in Reactome but not in background
    missing_genes = reactome_lib.unique_genes - background.genes
    
    print(f"\n🔍 Analysis Results:")
    print(f"   Reactome genes: {len(reactome_lib.unique_genes)}")
    print(f"   Background genes: {len(background.genes)}")
    print(f"   Genes in both: {len(reactome_lib.unique_genes & background.genes)}")
    print(f"   Missing genes: {len(missing_genes)}")
    
    if missing_genes:
        print(f"\n❌ Genes in Reactome but NOT in background ({len(missing_genes)} genes):")
        print("   " + ", ".join(sorted(list(missing_genes))[:20]))
        if len(missing_genes) > 20:
            print(f"   ... and {len(missing_genes) - 20} more")
        
        # Analyze which terms contain these missing genes
        print(f"\n📋 Terms containing missing genes:")
        missing_gene_terms = {}
        for term in reactome_lib.library:
            term_missing = set(term["genes"]) & missing_genes
            if term_missing:
                missing_gene_terms[term["name"]] = {
                    "missing_genes": term_missing,
                    "total_genes": len(term["genes"]),
                    "missing_count": len(term_missing)
                }
        
        # Sort by number of missing genes (descending)
        sorted_terms = sorted(missing_gene_terms.items(), 
                            key=lambda x: x[1]["missing_count"], reverse=True)
        
        print(f"   Terms with missing genes: {len(missing_gene_terms)}")
        print(f"\n🏆 Top 10 terms with most missing genes:")
        for i, (term_name, info) in enumerate(sorted_terms[:10], 1):
            print(f"   {i}. {term_name}")
            print(f"      Missing: {info['missing_count']}/{info['total_genes']} genes")
            print(f"      Missing genes: {', '.join(sorted(info['missing_genes']))}")
            print()
        
        # Check if there are any genes in background but not in Reactome
        extra_background_genes = background.genes - reactome_lib.unique_genes
        print(f"\n📈 Additional Analysis:")
        print(f"   Genes in background but NOT in Reactome: {len(extra_background_genes)}")
        
        # Show some examples of missing genes
        print(f"\n🔬 Sample Missing Genes Analysis:")
        sample_missing = list(missing_genes)[:10]
        for gene in sample_missing:
            # Find which terms contain this gene
            containing_terms = [term["name"] for term in reactome_lib.library if gene in term["genes"]]
            print(f"   {gene}: appears in {len(containing_terms)} terms")
            if containing_terms:
                print(f"      Example terms: {', '.join(containing_terms[:3])}")
                if len(containing_terms) > 3:
                    print(f"      ... and {len(containing_terms) - 3} more")
            print()
    
    else:
        print("\n✅ All Reactome genes are present in the background!")

def analyze_gene_overlap():
    """Analyze the overlap between Reactome and background genes."""
    
    print("\n\n📊 GENE OVERLAP ANALYSIS")
    print("=" * 50)
    
    # Load data
    reactome_lib = GeneSetLibrary(
        "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt",
        name="C2: CP: Reactome Pathways"
    )
    background = BackgroundGeneSet("data/backgrounds/all_genes.txt", name="All Genes")
    
    # Calculate overlaps
    reactome_genes = reactome_lib.unique_genes
    background_genes = background.genes
    
    intersection = reactome_genes & background_genes
    reactome_only = reactome_genes - background_genes
    background_only = background_genes - reactome_genes
    
    print(f"📈 Overlap Statistics:")
    print(f"   Reactome genes: {len(reactome_genes)}")
    print(f"   Background genes: {len(background_genes)}")
    print(f"   Intersection: {len(intersection)} ({len(intersection)/len(reactome_genes)*100:.1f}% of Reactome)")
    print(f"   Reactome only: {len(reactome_only)} ({len(reactome_only)/len(reactome_genes)*100:.1f}% of Reactome)")
    print(f"   Background only: {len(background_only)} ({len(background_only)/len(background_genes)*100:.1f}% of background)")
    
    # Calculate library-specific background size (as used in enrichment)
    library_background_size = len(intersection)
    print(f"\n🎯 Library-specific background size (for enrichment): {library_background_size}")

if __name__ == "__main__":
    find_missing_genes()
    analyze_gene_overlap()
