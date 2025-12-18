#!/usr/bin/env python3
"""
Example script showing how to use the p-value corrector with iGEA.

This demonstrates both approaches:
1. Post-processing (easier, recommended initially)
2. Integrated (more complex, better long-term)
"""

import sys
from pathlib import Path

# Add code directory to path
sys.path.insert(0, str(Path(__file__).parent / "code"))

from pvalue_corrector import load_corrector_from_permutations
from iter_enrichment import IterativeEnrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from background_gene_set import BackgroundGeneSet


def example_post_processing():
    """
    Approach 1: Post-processing correction
    
    Run iGEA normally, then apply correction afterward.
    """
    print("=" * 60)
    print("Approach 1: Post-Processing Correction")
    print("=" * 60)
    
    # 1. Initialize corrector (one-time setup)
    project_root = Path(__file__).parent
    permutation_file = project_root / "permutation_results" / "merged_permutation_results.tsv"
    
    if not permutation_file.exists():
        print(f"Error: Permutation file not found at {permutation_file}")
        print("Skipping correction example...")
        return
    
    print("\n1. Loading p-value corrector...")
    corrector = load_corrector_from_permutations(str(permutation_file))
    print(f"   ✓ Loaded corrector with {len([c for c in corrector.stratum_cdfs.values() if c is not None])} stratum CDFs")
    
    # 2. Setup iGEA (example - you'd use your actual data)
    print("\n2. Setting up iGEA...")
    # NOTE: Replace with your actual data loading
    # gene_set = GeneSet.load("your_genes.txt")
    # library = GeneSetLibrary.load("path/to/library.gmt")
    # background = BackgroundGeneSet.load("path/to/background.txt")
    
    print("   (Skipping actual iGEA run - replace with your data)")
    print("   Example code:")
    print("""
    iter_enrich = IterativeEnrichment(
        gene_set=gene_set,
        gene_set_library=library,
        background_gene_set=background,
        p_threshold=0.01,
        max_iterations=10,
    )
    """)
    
    # 3. Apply correction (example)
    print("\n3. Applying correction to results...")
    print("   Example code:")
    print("""
    # Convert to DataFrame
    df = iter_enrich.to_dataframe()
    
    # Add gene list size (needed for correction)
    df['gene_list_size'] = len(gene_set.genes)
    
    # Apply correction
    df_corrected = corrector.correct_dataframe(df)
    
    # Now df_corrected has 'corrected p-value' column
    print(df_corrected[['Term', 'iteration p-value', 'corrected p-value']].head())
    """)
    
    # 4. Demonstrate correction on a single p-value
    print("\n4. Example: Correcting a single p-value...")
    p_raw = 0.01
    library = "GO BP"
    iteration = 2
    term_size = 150
    overlap_size = 4
    gene_list_size = 73  # User has 73 genes (not in permutation data)
    
    p_corrected = corrector.correct_pvalue(
        p_raw, library, iteration, term_size, overlap_size, gene_list_size
    )
    
    print(f"   Raw p-value: {p_raw}")
    print(f"   Gene list size: {gene_list_size} → maps to {corrector._bin_gene_list_size(gene_list_size)}")
    print(f"   Corrected p-value: {p_corrected:.6f}")
    print(f"   Correction factor: {p_corrected / p_raw:.3f}x")
    
    print("\n✓ Post-processing approach complete!")


def example_integrated():
    """
    Approach 2: Integrated correction
    
    Pass corrector to IterativeEnrichment and correct p-values during computation.
    NOTE: This requires modifying IterativeEnrichment class (see integration guide).
    """
    print("\n" + "=" * 60)
    print("Approach 2: Integrated Correction")
    print("=" * 60)
    
    print("\n⚠️  This approach requires modifying IterativeEnrichment class.")
    print("   See PVALUE_CORRECTOR_INTEGRATION.md for details.")
    
    print("\nExample code (after modifying IterativeEnrichment):")
    print("""
    # 1. Initialize corrector
    corrector = load_corrector_from_permutations("permutation_results/merged_permutation_results.tsv")
    
    # 2. Run iGEA with corrector
    iter_enrich = IterativeEnrichment(
        gene_set=gene_set,
        gene_set_library=library,
        background_gene_set=background,
        pvalue_corrector=corrector,  # Pass corrector
        # ... other params
    )
    
    # Results now have corrected p-values
    df = iter_enrich.to_dataframe()
    # df['iteration p-value'] is already corrected
    """)
    
    print("\n✓ See PVALUE_CORRECTOR_INTEGRATION.md for full implementation details.")


def example_singleton():
    """
    Example: Using a singleton pattern for the corrector.
    """
    print("\n" + "=" * 60)
    print("Example: Singleton Pattern")
    print("=" * 60)
    
    print("""
# Create code/pvalue_corrector_singleton.py:

from code.pvalue_corrector import load_corrector_from_permutations
from pathlib import Path

_corrector = None

def get_corrector():
    \"\"\"Get or initialize the global corrector.\"\"\"
    global _corrector
    if _corrector is None:
        project_root = Path(__file__).parent.parent
        permutation_file = project_root / "permutation_results" / "merged_permutation_results.tsv"
        _corrector = load_corrector_from_permutations(str(permutation_file))
    return _corrector

# Usage:
from code.pvalue_corrector_singleton import get_corrector

corrector = get_corrector()  # Initialized once, reused everywhere
    """)


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("P-Value Corrector Usage Examples")
    print("=" * 60)
    
    # Run examples
    example_post_processing()
    example_integrated()
    example_singleton()
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("""
1. Initialize corrector once (takes ~10-30 seconds)
2. Use Approach 1 (post-processing) for easiest integration
3. Use Approach 2 (integrated) if you need corrected p-values during iteration
4. See PVALUE_CORRECTOR_INTEGRATION.md for detailed guide
    """)
