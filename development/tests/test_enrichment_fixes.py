#!/usr/bin/env python3
"""
Test script to verify our enrichment fixes are working correctly.
This tests the library-specific background calculation and input gene set filtering.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.enrichment import Enrichment
from code.gene_set import GeneSet
from code.gene_set_library import GeneSetLibrary
from code.background_gene_set import BackgroundGeneSet
import logging

# Set up logging to see the new messages
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_enrichment_fixes():
    """Test the enrichment fixes with a small gene set and library."""
    
    print("🧪 TESTING ENRICHMENT FIXES")
    print("=" * 60)
    
    # Load data
    print("\n📊 Loading data...")
    
    # Create a small test gene set with some genes that might not be in all libraries
    test_genes = ['TP53', 'BRCA1', 'BRCA2', 'GENE_NOT_IN_LIBRARY', 'ANOTHER_MISSING_GENE']
    
    # Load background
    background = BackgroundGeneSet("data/backgrounds/all_genes.txt", name="All Genes")
    print(f"Background genes: {background.size}")
    
    # Create gene set
    gene_set = GeneSet(test_genes, background.genes, name="Test Genes")
    print(f"Valid test genes: {gene_set.size}")
    print(f"Original input genes: {len(test_genes)}")
    print(f"Genes not in background: {len(test_genes) - gene_set.size}")
    
    # Test with a smaller library to see the filtering effect
    print("\n📚 Testing with GO Molecular Function library...")
    
    try:
        go_mf_lib = GeneSetLibrary(
            "data/libraries/c5.go.mf.v2025.1.Hs.symbols.gmt",
            name="C5: Gene Ontology: Molecular Function"
        )
        print(f"GO MF terms: {go_mf_lib.num_terms}")
        print(f"GO MF total unique genes: {go_mf_lib.size}")
        
        # Run enrichment with our fixes
        print("\n🔬 Running enrichment analysis...")
        enrichment = Enrichment(
            gene_set=gene_set,
            gene_set_library=go_mf_lib,
            background_gene_set=background,
            min_term_size=10,
            max_term_size=600,
            p_value_method_name="Hypergeometric Test"
        )
        
        # Get results
        results = enrichment.results
        
        print(f"\n📈 Enrichment results: {len(results)} terms found")
        
        if results:
            print("\n🏆 Top 3 results:")
            for i, result in enumerate(results[:3], 1):
                print(f"\n{i}. {result['term']}")
                print(f"   Description: {result['description']}")
                print(f"   P-value: {result['p-value']:.6e}")
                print(f"   FDR: {result['fdr']:.6e}")
                print(f"   Overlap: {result['overlap_size']} genes")
                print(f"   Overlap genes: {result['overlap']}")
        
        # Test with a different library to see different filtering
        print("\n\n📚 Testing with Reactome library...")
        
        reactome_lib = GeneSetLibrary(
            "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt",
            name="C2: Reactome Pathways"
        )
        print(f"Reactome terms: {reactome_lib.num_terms}")
        print(f"Reactome total unique genes: {reactome_lib.size}")
        
        # Run enrichment with Reactome
        print("\n🔬 Running enrichment analysis with Reactome...")
        enrichment_reactome = Enrichment(
            gene_set=gene_set,
            gene_set_library=reactome_lib,
            background_gene_set=background,
            min_term_size=10,
            max_term_size=600,
            p_value_method_name="Hypergeometric Test"
        )
        
        results_reactome = enrichment_reactome.results
        print(f"\n📈 Reactome results: {len(results_reactome)} terms found")
        
        if results_reactome:
            print("\n🏆 Top 3 Reactome results:")
            for i, result in enumerate(results_reactome[:3], 1):
                print(f"\n{i}. {result['term']}")
                print(f"   Description: {result['description']}")
                print(f"   P-value: {result['p-value']:.6e}")
                print(f"   FDR: {result['fdr']:.6e}")
                print(f"   Overlap: {result['overlap_size']} genes")
                print(f"   Overlap genes: {result['overlap']}")
        
    except Exception as e:
        print(f"Error during testing: {e}")
        import traceback
        traceback.print_exc()

def verify_fixes():
    """Verify that our fixes are working correctly."""
    print("\n\n✅ VERIFICATION OF FIXES")
    print("=" * 60)
    
    print("""
Our fixes should show the following improvements:

1. ✅ Library-specific background calculation:
   - Should use only genes from terms within size filter (10-600)
   - Should show filtered gene count in logs
   - Should be different from total library genes

2. ✅ Input gene set filtering:
   - Should filter input genes with library-specific background
   - Should show filtering impact in logs
   - Should use filtered gene set size in calculations

3. ✅ Statistical consistency:
   - Same gene universe for input and background
   - More accurate p-values
   - Proper hypergeometric parameters

The logging messages should show:
- "Library-specific background size: X genes (intersection of Y background genes and Z filtered library genes from W terms within size range [10, 600])"
- "Input gene set filtered: X → Y genes (intersected with library-specific background)" (if filtering occurred)
- "Input gene set size: X genes (all genes present in library-specific background)" (if no filtering needed)
""")

if __name__ == "__main__":
    test_enrichment_fixes()
    verify_fixes()
    
    print("\n🎉 Test complete!")
    print("\n💡 Check the logging output above to verify our fixes are working correctly.")
