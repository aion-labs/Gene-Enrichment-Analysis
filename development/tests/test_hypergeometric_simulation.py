#!/usr/bin/env python3
"""
Test script to verify hypergeometric p-value calculation for Reactome with hypoxia gene list.
This shows the exact values passed to scipy.stats.hypergeom.sf function.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.gene_set import GeneSet
from code.gene_set_library import GeneSetLibrary
from code.background_gene_set import BackgroundGeneSet
from scipy.stats import hypergeom
import json

def simulate_hypergeometric_test():
    """Simulate the hypergeometric test calculation."""
    
    print("🔬 HYPERGEOMETRIC TEST SIMULATION")
    print("=" * 60)
    
    # Load data
    print("\n📊 Loading data...")
    
    # Load hypoxia genes
    try:
        with open('data/gene_lists/hypoxia-genes.symbols.txt', 'r') as f:
            hypoxia_genes = [line.strip() for line in f if line.strip()]
        print(f"Hypoxia genes: {len(hypoxia_genes)}")
    except FileNotFoundError:
        print("Error: hypoxia-genes.symbols.txt not found")
        return
    
    # Create background gene set from file
    try:
        background = BackgroundGeneSet("data/backgrounds/all_genes.txt", name="All Genes")
        print(f"Background genes: {background.size}")
    except FileNotFoundError:
        print("Error: all_genes.txt not found")
        return
    
    # Create gene set
    gene_set = GeneSet(hypoxia_genes, background.genes, name="Hypoxia Genes")
    print(f"Valid hypoxia genes: {gene_set.size}")
    
    # Load Reactome library
    print("\n📚 Loading Reactome library...")
    try:
        reactome_lib = GeneSetLibrary(
            "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt",
            name="C2: CP: Reactome Pathways"
        )
        print(f"Reactome terms: {reactome_lib.num_terms}")
        print(f"Reactome unique genes: {reactome_lib.size}")
    except Exception as e:
        print(f"Error loading Reactome library: {e}")
        return
    
    # Calculate library-specific background size
    library_background_size = len(background.genes & reactome_lib.unique_genes)
    print(f"\n🎯 Library-specific background size: {library_background_size}")
    print(f"(Intersection of {background.size} background genes and {reactome_lib.size} Reactome genes)")
    
    # Test a few terms to show the calculation
    print("\n🧮 Testing hypergeometric calculation for sample terms:")
    print("-" * 60)
    
    # Filter terms by size (similar to the actual implementation)
    min_term_size = 10
    max_term_size = 600
    
    # Find terms that have overlap with hypoxia genes
    terms_with_overlap = []
    for term in reactome_lib.library:
        if min_term_size <= term["size"] <= max_term_size:
            term_genes = set(term["genes"])
            overlap = gene_set.genes & term_genes
            if len(overlap) > 0:  # Only include terms with overlap
                terms_with_overlap.append((term, len(overlap)))
    
    # Sort by overlap size (descending) and take top 5
    terms_with_overlap.sort(key=lambda x: x[1], reverse=True)
    test_terms = [term for term, _ in terms_with_overlap[:5]]
    
    if not test_terms:
        print("No terms found with overlap. Testing with all terms...")
        test_terms = [term for term in reactome_lib.library 
                      if min_term_size <= term["size"] <= max_term_size][:5]
    
    for i, term in enumerate(test_terms, 1):
        print(f"\n📋 Term {i}: {term['name']}")
        print(f"   Description: {term['description']}")
        
        # Calculate overlap
        term_genes = set(term["genes"])
        overlap = gene_set.genes & term_genes
        n_overlap = len(overlap)
        n_term_genes = len(term_genes)
        
        print(f"   Term size: {n_term_genes}")
        print(f"   Overlap: {n_overlap}")
        print(f"   Overlap genes: {sorted(list(overlap))}")
        
        # Show hypergeometric parameters
        print(f"\n   📐 Hypergeometric parameters:")
        print(f"      M (total population): {library_background_size}")
        print(f"      n (successes in population): {n_term_genes}")
        print(f"      N (sample size): {gene_set.size}")
        print(f"      k (successes in sample): {n_overlap}")
        
        # Calculate p-value using hypergeom.sf
        p_value = hypergeom.sf(
            n_overlap - 1,  # k - 1 (survival function gives P(X > k-1) = P(X >= k))
            library_background_size,  # M
            n_term_genes,  # n
            gene_set.size  # N
        )
        
        print(f"\n   🧮 hypergeom.sf({n_overlap - 1}, {library_background_size}, {n_term_genes}, {gene_set.size})")
        print(f"   P-value: {p_value:.6e}")
        
        # Also show the exact scipy call
        print(f"   📞 Exact scipy call: hypergeom.sf(k={n_overlap-1}, M={library_background_size}, n={n_term_genes}, N={gene_set.size})")
        
        # Show what this means
        print(f"   📖 Interpretation: P(X >= {n_overlap}) where X ~ Hypergeometric(M={library_background_size}, n={n_term_genes}, N={gene_set.size})")
        
        if i < len(test_terms):
            print("   " + "-" * 50)

def verify_hypergeometric_formula():
    """Verify the hypergeometric formula and parameters."""
    print("\n\n🔍 HYPERGEOMETRIC FORMULA VERIFICATION")
    print("=" * 60)
    
    print("""
The hypergeometric distribution models the probability of drawing k successes 
from a population of size M containing n successes, when drawing N items.

Parameters:
- M: Total population size (library_background_size)
- n: Number of successes in population (term size)
- N: Sample size (gene set size)
- k: Number of successes in sample (overlap size)

Formula: P(X = k) = C(n,k) * C(M-n, N-k) / C(M,N)

For enrichment analysis, we want P(X >= k), which is calculated using the 
survival function: hypergeom.sf(k-1, M, n, N)

This gives us the probability of observing k or more successes by chance.
""")

if __name__ == "__main__":
    simulate_hypergeometric_test()
    verify_hypergeometric_formula()
    
    print("\n✅ Simulation complete!")
    print("\n💡 Key points to verify:")
    print("1. M (library_background_size) should be the intersection of background and library genes")
    print("2. n (term size) should be the number of genes in the specific term")
    print("3. N (gene set size) should be the number of genes in your input set")
    print("4. k (overlap) should be the number of genes that appear in both term and input set")
    print("5. The survival function gives P(X >= k) = 1 - P(X < k)")
