#!/usr/bin/env python3
"""
Simple test to verify ORA vs iGEA comparison statistics function.
Tests the core logic directly with mock data structures.
"""

import sys
from pathlib import Path

# Add code directory to path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))

from ui.rendering import calculate_ora_igea_comparison


class MockEnrichment:
    """Mock Enrichment object for testing."""
    def __init__(self, results):
        self.results = results


class MockIterativeEnrichment:
    """Mock IterativeEnrichment object for testing."""
    def __init__(self, results):
        self.results = results


def test_basic_comparison():
    """Test basic comparison with mock data."""
    print("=" * 60)
    print("Testing calculate_ora_igea_comparison with mock data")
    print("=" * 60)
    
    # Create mock regular enrichment results
    regular_enrichments = {
        "Test Library": MockEnrichment([
            {
                "term": "TERM1",
                "overlap": ["GENE1", "GENE2", "GENE3"],
                "p-value": 0.001,  # Significant
                "fdr": 0.01,
            },
            {
                "term": "TERM2",
                "overlap": ["GENE3", "GENE4", "GENE5"],
                "p-value": 0.005,  # Significant
                "fdr": 0.02,
            },
            {
                "term": "TERM3",
                "overlap": ["GENE6", "GENE7"],
                "p-value": 0.1,  # Not significant
                "fdr": 0.5,
            }
        ])
    }
    
    # Create mock iterative enrichment results
    iterative_enrichments = {
        "Test Library": MockIterativeEnrichment([
            {
                "Iteration": 1,
                "Term": "TERM1",
                "iteration p-value": 0.001,  # Significant
                "Genes removed for next iteration": ["GENE1", "GENE2", "GENE3"],
            },
            {
                "Iteration": 2,
                "Term": "TERM2",
                "iteration p-value": 0.005,  # Significant
                "Genes removed for next iteration": ["GENE4", "GENE5"],
            },
            {
                "Iteration": 3,
                "Term": "TERM4",
                "iteration p-value": 0.002,  # Significant
                "Genes removed for next iteration": ["GENE8", "GENE9"],
            }
        ])
    }
    
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
    print(f"  Genes - ORA Only: {stats['genes']['ora_only']}")
    print(f"  Genes - iGEA Only: {stats['genes']['igea_only']}")
    print(f"  Terms - ORA Total: {stats['terms']['ora_total']}")
    print(f"  Terms - iGEA Total: {stats['terms']['igea_total']}")
    print(f"  Terms - Overlap: {stats['terms']['overlap']}")
    print(f"  Terms - ORA Only: {stats['terms']['ora_only']}")
    print(f"  Terms - iGEA Only: {stats['terms']['igea_only']}")
    
    # Verify expected results
    assert stats['genes']['ora_total'] == 5, f"Expected 5 ORA genes, got {stats['genes']['ora_total']}"
    assert stats['genes']['igea_total'] == 7, f"Expected 7 iGEA genes, got {stats['genes']['igea_total']}"
    assert stats['genes']['overlap'] == 5, f"Expected 5 overlapping genes, got {stats['genes']['overlap']}"
    assert stats['terms']['ora_total'] == 2, f"Expected 2 ORA terms, got {stats['terms']['ora_total']}"
    assert stats['terms']['igea_total'] == 3, f"Expected 3 iGEA terms, got {stats['terms']['igea_total']}"
    assert stats['terms']['overlap'] == 2, f"Expected 2 overlapping terms, got {stats['terms']['overlap']}"
    
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


def test_condition_check():
    """Test the condition that should be checked in streamlit_app.py"""
    print("\n" + "=" * 60)
    print("Testing condition check (as in streamlit_app.py)")
    print("=" * 60)
    
    # Simulate state.enrich and state.iter_enrich
    state_enrich = {
        "Library1": MockEnrichment([{"term": "TERM1", "overlap": ["GENE1"], "p-value": 0.001, "fdr": 0.01}])
    }
    state_iter_enrich = {
        "Library1": MockIterativeEnrichment([{"Term": "TERM1", "iteration p-value": 0.001, "Genes removed for next iteration": ["GENE1"]}])
    }
    
    # This is the condition from streamlit_app.py
    condition_met = state_enrich and state_iter_enrich and len(state_enrich) > 0 and len(state_iter_enrich) > 0
    
    print(f"  state.enrich exists: {bool(state_enrich)}")
    print(f"  state.iter_enrich exists: {bool(state_iter_enrich)}")
    print(f"  len(state.enrich) > 0: {len(state_enrich) > 0}")
    print(f"  len(state.iter_enrich) > 0: {len(state_iter_enrich) > 0}")
    print(f"  Condition met: {condition_met}")
    
    if condition_met:
        comparison = calculate_ora_igea_comparison(
            state_enrich,
            state_iter_enrich,
            p_threshold=0.01,
            fdr_threshold=0.05
        )
        print(f"  ✅ Comparison calculated successfully!")
        print(f"  Genes ORA: {comparison['stats']['genes']['ora_total']}")
        print(f"  Genes iGEA: {comparison['stats']['genes']['igea_total']}")
    else:
        print("  ❌ Condition not met - comparison would not be displayed")
    
    assert condition_met, "Condition should be met"
    print("✅ Condition check test passed!")


def main():
    """Run all tests."""
    print("\n" + "🧪 " * 20)
    print("ORA vs iGEA Comparison Statistics Test Suite (Simple)")
    print("🧪 " * 20 + "\n")
    
    try:
        # Test 1: Basic comparison
        comparison = test_basic_comparison()
        
        # Test 2: Empty results
        test_empty_results()
        
        # Test 3: Condition check
        test_condition_check()
        
        print("\n" + "=" * 60)
        print("✅ ALL TESTS PASSED!")
        print("=" * 60)
        print("\nThe comparison statistics function is working correctly.")
        print("\nTo see it in Streamlit:")
        print("  1. Run Regular enrichment first")
        print("  2. Then run Iterative enrichment")
        print("  3. The comparison should appear at the bottom of iterative results")
        print("  4. Check that both state.enrich and state.iter_enrich are populated")
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

