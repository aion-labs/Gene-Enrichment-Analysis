#!/usr/bin/env python3
"""
Test library name matching for benchmarking in CLI and Streamlit.

This script verifies that all 11 libraries with permutation data
can be correctly matched between user library names and parquet library names.
"""

import sys
from pathlib import Path

# Add code directory to path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from parallel_null_distribution import (
    find_intersection_libraries,
    normalize_library_name,
    get_available_libraries_from_parquet
)

# All 11 libraries that should be in parquet files
# These are the user-friendly names that would be used in iGEA
USER_LIBRARY_NAMES = [
    "Reactome",
    "KEGG",
    "GO BP",
    "GO MF",
    "GO CC",
    "BioCarta",
    "Canonical pathways",
    "KEGG Medicus",
    "Pathway Interaction Database",
    "WikiPathways",
    "Hallmark",
]

# Alternative names that might be used (from Streamlit/CLI)
ALTERNATIVE_NAMES = {
    "Reactome": ["C2: CP: Reactome Pathways", "C2: Reactome Pathways"],
    "KEGG": ["C2: CP: KEGG Legacy", "C2: KEGG Legacy"],
    "GO BP": ["C5: GO: Biological Process", "C5: Gene Ontology: Biological Process"],
    "GO MF": ["C5: GO: Molecular Function", "C5: Gene Ontology: Molecular Function"],
    "GO CC": ["C5: GO: Cellular Component", "C5: Gene Ontology: Cellular Component"],
    "BioCarta": ["C2: CP: BioCarta"],
    "Canonical pathways": ["C2: CP: Canonical pathways", "C2: Canonical pathways"],
    "KEGG Medicus": ["C2: CP: KEGG MEDICUS", "C2: KEGG MEDICUS"],
    "Pathway Interaction Database": ["C2: CP: Pathway Interaction Database", "C2: Pathway Interaction Database"],
    "WikiPathways": ["C2: CP: WikiPathways", "C2: WikiPathways"],
    "Hallmark": ["H: Hallmark Gene Sets"],
}


def test_normalization():
    """Test that normalization works correctly for all libraries."""
    print("=" * 80)
    print("TEST 1: Library Name Normalization")
    print("=" * 80)
    print()
    
    parquet_dir = PROJECT_ROOT / "permutations" / "permutation_cluster_statistics_parquet"
    available_libraries = get_available_libraries_from_parquet(parquet_dir)
    
    print(f"Available libraries in parquet: {len(available_libraries)}")
    print(f"Libraries: {sorted(available_libraries)}")
    print()
    
    all_passed = True
    
    for user_name in USER_LIBRARY_NAMES:
        norm_user = normalize_library_name(user_name)
        print(f"User library: '{user_name}'")
        print(f"  Normalized: '{norm_user}'")
        
        # Try to match with parquet libraries
        matched = False
        for parquet_lib in available_libraries:
            norm_parquet = normalize_library_name(parquet_lib)
            if (norm_user == norm_parquet or 
                norm_user in norm_parquet.lower() or 
                norm_parquet in norm_user.lower()):
                print(f"  ✓ Matches: '{parquet_lib}' (normalized: '{norm_parquet}')")
                matched = True
                break
        
        if not matched:
            print(f"  ✗ NO MATCH FOUND")
            all_passed = False
        print()
    
    return all_passed


def test_find_intersection():
    """Test find_intersection_libraries function (used in CLI and Streamlit)."""
    print("=" * 80)
    print("TEST 2: find_intersection_libraries (CLI/Streamlit matching)")
    print("=" * 80)
    print()
    
    parquet_dir = PROJECT_ROOT / "permutations" / "permutation_cluster_statistics_parquet"
    
    all_passed = True
    
    # Test with user library names
    print("Testing with user library names:")
    libraries_with_data, libraries_without_data = find_intersection_libraries(
        USER_LIBRARY_NAMES,
        parquet_dir,
        use_all_available=False
    )
    
    print(f"  Libraries with data: {len(libraries_with_data)}")
    for lib in libraries_with_data:
        print(f"    - {lib}")
    
    print(f"  Libraries without data: {len(libraries_without_data)}")
    for lib in libraries_without_data:
        print(f"    - {lib}")
    
    if len(libraries_with_data) < 11:
        print(f"  ✗ WARNING: Only {len(libraries_with_data)}/11 libraries matched!")
        all_passed = False
    else:
        print(f"  ✓ All 11 libraries matched!")
    
    print()
    
    # Test with alternative names (as might appear in Streamlit/CLI)
    print("Testing with alternative library names (as might appear in Streamlit/CLI):")
    alternative_user_names = []
    for base_name in USER_LIBRARY_NAMES:
        # Use alternative names if available, otherwise use base name
        if base_name in ALTERNATIVE_NAMES:
            alternative_user_names.extend(ALTERNATIVE_NAMES[base_name])
        else:
            alternative_user_names.append(base_name)
    
    # Remove duplicates while preserving order
    seen = set()
    unique_alt_names = []
    for name in alternative_user_names:
        if name not in seen:
            seen.add(name)
            unique_alt_names.append(name)
    
    libraries_with_data_alt, libraries_without_data_alt = find_intersection_libraries(
        unique_alt_names,
        parquet_dir,
        use_all_available=False
    )
    
    print(f"  Libraries with data: {len(libraries_with_data_alt)}")
    for lib in libraries_with_data_alt:
        print(f"    - {lib}")
    
    print(f"  Libraries without data: {len(libraries_without_data_alt)}")
    for lib in libraries_without_data_alt:
        print(f"    - {lib}")
    
    if len(libraries_with_data_alt) < 11:
        print(f"  ✗ WARNING: Only {len(libraries_with_data_alt)}/11 libraries matched with alternative names!")
        all_passed = False
    else:
        print(f"  ✓ All 11 libraries matched with alternative names!")
    
    print()
    
    return all_passed


def test_mapping_logic():
    """Test the mapping logic used in CLI and Streamlit benchmarking."""
    print("=" * 80)
    print("TEST 3: User-to-Parquet Mapping Logic (as used in benchmarking)")
    print("=" * 80)
    print()
    
    parquet_dir = PROJECT_ROOT / "permutations" / "permutation_cluster_statistics_parquet"
    available_libraries = get_available_libraries_from_parquet(parquet_dir)
    
    # Simulate the mapping logic from CLI/Streamlit
    user_to_parquet_mapping = {}
    
    for user_lib in USER_LIBRARY_NAMES:
        normalized_user = normalize_library_name(user_lib)
        for parquet_lib in available_libraries:
            normalized_parquet = normalize_library_name(parquet_lib)
            if (normalized_user == normalized_parquet or 
                normalized_user in normalized_parquet.lower() or 
                normalized_parquet in normalized_user.lower()):
                user_to_parquet_mapping[user_lib] = parquet_lib
                break
    
    print(f"Mapping results: {len(user_to_parquet_mapping)}/{len(USER_LIBRARY_NAMES)} libraries matched")
    print()
    
    all_passed = True
    for user_lib, parquet_lib in sorted(user_to_parquet_mapping.items()):
        print(f"  ✓ '{user_lib}' → '{parquet_lib}'")
    
    print()
    
    missing = set(USER_LIBRARY_NAMES) - set(user_to_parquet_mapping.keys())
    if missing:
        print(f"  ✗ Missing mappings: {missing}")
        all_passed = False
    else:
        print(f"  ✓ All {len(USER_LIBRARY_NAMES)} libraries have mappings!")
    
    print()
    
    return all_passed


def main():
    """Run all tests."""
    print()
    print("=" * 80)
    print("LIBRARY NAME MATCHING TEST")
    print("Testing that all 11 libraries match correctly for benchmarking")
    print("=" * 80)
    print()
    
    results = []
    
    # Test 1: Normalization
    results.append(("Normalization", test_normalization()))
    
    # Test 2: find_intersection_libraries
    results.append(("find_intersection_libraries", test_find_intersection()))
    
    # Test 3: Mapping logic
    results.append(("Mapping Logic", test_mapping_logic()))
    
    # Summary
    print("=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    print()
    
    all_passed = True
    for test_name, passed in results:
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"  {test_name}: {status}")
        if not passed:
            all_passed = False
    
    print()
    if all_passed:
        print("✓ ALL TESTS PASSED - Library matching works correctly!")
    else:
        print("✗ SOME TESTS FAILED - Please check the output above")
    
    print()
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())

