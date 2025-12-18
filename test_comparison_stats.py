#!/usr/bin/env python3
"""
Test script to verify ORA vs iGEA comparison statistics implementation.
This script tests the comparison function directly without Streamlit.
"""

import sys
from pathlib import Path

# Add code directory to path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))

from enrichment import Enrichment
from iter_enrichment import IterativeEnrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from background_gene_set import BackgroundGeneSet
from ui.rendering import calculate_ora_igea_comparison, render_ora_igea_comparison

def create_mock_enrichment_results():
    """Create mock enrichment results for testing."""
    import tempfile
    from pathlib import Path
    
    # Create a temporary background file
    background_genes = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5", 
                        "GENE6", "GENE7", "GENE8", "GENE9", "GENE10"]
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        f.write('\n'.join(background_genes))
        temp_bg_file = f.name
    
    try:
        # Create a simple gene set
        genes = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
        background_genes_set = set(background_genes)
        
        gene_set = GeneSet(genes, background_genes_set, hgcn=False, format=False)
        background_set = BackgroundGeneSet(temp_bg_file, "test_background")
        
        # Create a simple library
        library_terms = [
            {
                "term": "TERM1",
                "description": "Test term 1",
                "genes": ["GENE1", "GENE2", "GENE3"],
                "size": 3
            },
            {
                "term": "TERM2",
                "description": "Test term 2",
                "genes": ["GENE3", "GENE4", "GENE5"],
                "size": 3
            }
        ]
        
        library = GeneSetLibrary("Test Library", library_terms)
        
        # Create mock regular enrichment results
        regular_enrich = Enrichment(
            gene_set=gene_set,
            gene_set_library=library,
            background_gene_set=background_set,
            min_term_size=1,
            max_term_size=10,
            p_value_method_name="Fisher's Exact Test"
        )
        
        # Manually set results for regular enrichment
        regular_enrich.results = [
            {
                "term": "TERM1",
                "description": "Test term 1",
                "overlap": ["GENE1", "GENE2", "GENE3"],
                "overlap_size": "3/3",
                "p-value": 0.001,  # Significant
                "fdr": 0.01,
                "rank": 1
            },
            {
                "term": "TERM2",
                "description": "Test term 2",
                "overlap": ["GENE3", "GENE4", "GENE5"],
                "overlap_size": "3/3",
                "p-value": 0.005,  # Significant
                "fdr": 0.02,
                "rank": 2
            }
        ]
        
        # Create mock iterative enrichment - we'll create a minimal one
        # Since IterativeEnrichment runs computation in __init__, we'll mock it differently
        class MockIterativeEnrichment:
            def __init__(self):
                self.results = [
                    {
                        "Iteration": 1,
                        "Term": "TERM1",
                        "Description": "Test term 1",
                        "Library": "Test Library",
                        "iteration p-value": 0.001,  # Significant
                        "iteration overlapping genes": "3/3",
                        "Genes removed for next iteration": ["GENE1", "GENE2", "GENE3"],
                    },
                    {
                        "Iteration": 2,
                        "Term": "TERM2",
                        "Description": "Test term 2",
                        "Library": "Test Library",
                        "iteration p-value": 0.005,  # Significant
                        "iteration overlapping genes": "2/3",
                        "Genes removed for next iteration": ["GENE4", "GENE5"],
                    }
                ]
        
        iter_enrich = MockIterativeEnrichment()
        
        return {
            "regular": {"Test Library": regular_enrich},
            "iterative": {"Test Library": iter_enrich}
        }
    finally:
        # Clean up temp file
        Path(temp_bg_file).unlink(missing_ok=True)


def test_calculate_comparison():
    """Test the calculate_ora_igea_comparison function."""
    print("=" * 60)
    print("Testing calculate_ora_igea_comparison function")
    print("=" * 60)
    
    mock_data = create_mock_enrichment_results()
    regular_enrichments = mock_data["regular"]
    iterative_enrichments = mock_data["iterative"]
    
    # Test the calculation
    comparison = calculate_ora_igea_comparison(
        regular_enrichments,
        iterative_enrichments,
        p_threshold=0.01,
        fdr_threshold=0.05
    )
    
    # Print results
    print("\n📊 Comparison Results:")
    print(f"  Genes ORA: {len(comparison['genes_ora'])} - {sorted(comparison['genes_ora'])}")
    print(f"  Genes iGEA: {len(comparison['genes_igea'])} - {sorted(comparison['genes_igea'])}")
    print(f"  Genes Overlap: {len(comparison['genes_overlap'])} - {sorted(comparison['genes_overlap'])}")
    print(f"  Genes ORA Only: {len(comparison['genes_ora_only'])} - {sorted(comparison['genes_ora_only'])}")
    print(f"  Genes iGEA Only: {len(comparison['genes_igea_only'])} - {sorted(comparison['genes_igea_only'])}")
    
    print(f"\n  Terms ORA: {len(comparison['terms_ora'])} - {sorted(comparison['terms_ora'])}")
    print(f"  Terms iGEA: {len(comparison['terms_igea'])} - {sorted(comparison['terms_igea'])}")
    print(f"  Terms Overlap: {len(comparison['terms_overlap'])} - {sorted(comparison['terms_overlap'])}")
    print(f"  Terms ORA Only: {len(comparison['terms_ora_only'])} - {sorted(comparison['terms_ora_only'])}")
    print(f"  Terms iGEA Only: {len(comparison['terms_igea_only'])} - {sorted(comparison['terms_igea_only'])}")
    
    print("\n📈 Statistics:")
    stats = comparison["stats"]
    print(f"  Genes - ORA Total: {stats['genes']['ora_total']}")
    print(f"  Genes - iGEA Total: {stats['genes']['igea_total']}")
    print(f"  Genes - Overlap: {stats['genes']['overlap']}")
    print(f"  Terms - ORA Total: {stats['terms']['ora_total']}")
    print(f"  Terms - iGEA Total: {stats['terms']['igea_total']}")
    print(f"  Terms - Overlap: {stats['terms']['overlap']}")
    
    # Verify expected results
    assert len(comparison['genes_ora']) > 0, "Should have genes from ORA"
    assert len(comparison['genes_igea']) > 0, "Should have genes from iGEA"
    assert len(comparison['terms_ora']) > 0, "Should have terms from ORA"
    assert len(comparison['terms_igea']) > 0, "Should have terms from iGEA"
    
    print("\n✅ All assertions passed!")
    return comparison


def test_empty_results():
    """Test with empty results."""
    print("\n" + "=" * 60)
    print("Testing with empty results")
    print("=" * 60)
    
    comparison = calculate_ora_igea_comparison(
        {},
        {},
        p_threshold=0.01,
        fdr_threshold=0.05
    )
    
    assert comparison['stats']['genes']['ora_total'] == 0
    assert comparison['stats']['genes']['igea_total'] == 0
    assert comparison['stats']['terms']['ora_total'] == 0
    assert comparison['stats']['terms']['igea_total'] == 0
    
    print("✅ Empty results test passed!")


def test_streamlit_integration():
    """Test that the function can be called in a Streamlit-like context."""
    print("\n" + "=" * 60)
    print("Testing Streamlit integration (mock)")
    print("=" * 60)
    
    # Mock streamlit module
    class MockStreamlit:
        def markdown(self, text):
            print(f"  [MARKDOWN] {text}")
        
        def metric(self, label, value):
            print(f"  [METRIC] {label}: {value}")
        
        def caption(self, text):
            print(f"  [CAPTION] {text}")
        
        def columns(self, n):
            return [self] * n
        
        def info(self, text):
            print(f"  [INFO] {text}")
        
        def __enter__(self):
            return self
        
        def __exit__(self, *args):
            pass
    
    # Temporarily replace streamlit module
    import ui.rendering as rendering_module
    original_st = rendering_module.st
    rendering_module.st = MockStreamlit()
    
    try:
        mock_data = create_mock_enrichment_results()
        print("\n  Calling render_ora_igea_comparison...")
        render_ora_igea_comparison(
            mock_data["regular"],
            mock_data["iterative"],
            p_threshold=0.01,
            fdr_threshold=0.05
        )
        print("  ✅ render_ora_igea_comparison executed successfully!")
    except Exception as e:
        print(f"  ❌ Error: {e}")
        raise
    finally:
        # Restore original streamlit
        rendering_module.st = original_st


def main():
    """Run all tests."""
    print("\n" + "🧪 " * 20)
    print("ORA vs iGEA Comparison Statistics Test Suite")
    print("🧪 " * 20 + "\n")
    
    try:
        # Test 1: Basic calculation
        comparison = test_calculate_comparison()
        
        # Test 2: Empty results
        test_empty_results()
        
        # Test 3: Streamlit integration (mock)
        test_streamlit_integration()
        
        print("\n" + "=" * 60)
        print("✅ ALL TESTS PASSED!")
        print("=" * 60)
        print("\nThe comparison statistics function is working correctly.")
        print("If you don't see it in Streamlit, check:")
        print("  1. Both state.enrich and state.iter_enrich are populated")
        print("  2. Both have results (not empty)")
        print("  3. The condition check in streamlit_app.py is being met")
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

