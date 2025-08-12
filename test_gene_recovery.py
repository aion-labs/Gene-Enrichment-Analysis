#!/usr/bin/env python3
"""
Script to test if we can recover the missing genes using gene validation methods.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.gene_set_library import GeneSetLibrary
from code.background_gene_set import BackgroundGeneSet
from code.gene_converter import GeneConverter

def test_gene_recovery():
    """Test if we can recover missing genes using gene validation methods."""
    
    print("🔍 TESTING GENE RECOVERY")
    print("=" * 50)
    
    # Load Reactome library and background
    print("📚 Loading data...")
    reactome_lib = GeneSetLibrary(
        "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt",
        name="C2: CP: Reactome Pathways"
    )
    background = BackgroundGeneSet("data/backgrounds/all_genes.txt", name="All Genes")
    
    # Find missing genes
    missing_genes = reactome_lib.unique_genes - background.genes
    print(f"Missing genes: {len(missing_genes)}")
    
    # Initialize gene converter
    converter = GeneConverter()
    
    print(f"\n🧬 Testing gene recovery for {len(missing_genes)} missing genes:")
    print("-" * 60)
    
    recovered_genes = {}
    unrecovered_genes = []
    
    for gene in sorted(missing_genes):
        print(f"\n🔍 Testing: {gene}")
        
        # Test direct validation
        direct_result = converter.validate_and_map_symbol(gene)
        
        # Test if it's an Entrez ID
        entrez_result = None
        if converter.is_entrez_id(gene):
            entrez_result = converter.get_symbol(gene)
        
        # Test synonyms and history
        synonyms_result = None
        history_result = None
        
        # Check if it's in synonyms
        if gene.upper() in converter.synonyms_to_symbol:
            synonyms_result = converter.synonyms_to_symbol[gene.upper()]
        
        # Check if it's in gene history
        if gene.upper() in converter.old_to_current_symbol:
            history_result = converter.old_to_current_symbol[gene.upper()]
        
        # Report results
        print(f"   Direct validation: {direct_result}")
        print(f"   Entrez ID check: {entrez_result}")
        print(f"   Synonyms mapping: {synonyms_result}")
        print(f"   History mapping: {history_result}")
        
        # Determine if we can recover this gene
        recovered_symbol = None
        recovery_method = None
        
        if direct_result:
            recovered_symbol = direct_result
            recovery_method = "direct_validation"
        elif entrez_result:
            recovered_symbol = entrez_result
            recovery_method = "entrez_id"
        elif synonyms_result:
            recovered_symbol = synonyms_result
            recovery_method = "synonyms"
        elif history_result:
            recovered_symbol = history_result
            recovery_method = "gene_history"
        
        if recovered_symbol:
            # Check if the recovered symbol is in the background
            if recovered_symbol in background.genes:
                recovered_genes[gene] = {
                    "recovered_symbol": recovered_symbol,
                    "method": recovery_method,
                    "in_background": True
                }
                print(f"   ✅ RECOVERED: {gene} -> {recovered_symbol} (via {recovery_method})")
            else:
                recovered_genes[gene] = {
                    "recovered_symbol": recovered_symbol,
                    "method": recovery_method,
                    "in_background": False
                }
                print(f"   ⚠️  PARTIALLY RECOVERED: {gene} -> {recovered_symbol} (via {recovery_method}) but not in background")
        else:
            unrecovered_genes.append(gene)
            print(f"   ❌ NOT RECOVERED: {gene}")
    
    # Summary
    print(f"\n\n📊 RECOVERY SUMMARY")
    print("=" * 50)
    print(f"Total missing genes: {len(missing_genes)}")
    print(f"Recovered genes: {len(recovered_genes)}")
    print(f"Unrecovered genes: {len(unrecovered_genes)}")
    
    if recovered_genes:
        print(f"\n✅ RECOVERED GENES ({len(recovered_genes)}):")
        for gene, info in recovered_genes.items():
            status = "✅" if info["in_background"] else "⚠️"
            print(f"   {status} {gene} -> {info['recovered_symbol']} (via {info['method']})")
    
    if unrecovered_genes:
        print(f"\n❌ UNRECOVERED GENES ({len(unrecovered_genes)}):")
        for gene in unrecovered_genes:
            print(f"   ❌ {gene}")
    
    # Analyze recovery by method
    if recovered_genes:
        print(f"\n📈 RECOVERY BY METHOD:")
        method_counts = {}
        for info in recovered_genes.values():
            method = info["method"]
            method_counts[method] = method_counts.get(method, 0) + 1
        
        for method, count in method_counts.items():
            print(f"   {method}: {count} genes")

def test_mitochondrial_genes():
    """Specifically test mitochondrial genes which are likely to be problematic."""
    
    print(f"\n\n🔬 MITOCHONDRIAL GENE ANALYSIS")
    print("=" * 50)
    
    converter = GeneConverter()
    background = BackgroundGeneSet("data/backgrounds/all_genes.txt", name="All Genes")
    
    # Mitochondrial genes from the missing list
    mt_genes = [
        "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO3", "MT-CYB",
        "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6",
        "MT-RNR1", "MT-RNR2"
    ]
    
    print(f"Testing {len(mt_genes)} mitochondrial genes:")
    
    for gene in mt_genes:
        print(f"\n🔍 {gene}:")
        
        # Check if it's in the background directly
        in_background = gene in background.genes
        print(f"   In background: {in_background}")
        
        # Try validation
        validated = converter.validate_and_map_symbol(gene)
        print(f"   Validated: {validated}")
        
        # Check if it's in synonyms
        in_synonyms = gene.upper() in converter.synonyms_to_symbol
        if in_synonyms:
            synonym_result = converter.synonyms_to_symbol[gene.upper()]
            print(f"   In synonyms: {synonym_result}")
        
        # Check if it's in gene history
        in_history = gene.upper() in converter.old_to_current_symbol
        if in_history:
            history_result = converter.old_to_current_symbol[gene.upper()]
            print(f"   In history: {history_result}")

if __name__ == "__main__":
    test_gene_recovery()
    test_mitochondrial_genes()
