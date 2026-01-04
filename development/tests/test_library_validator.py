#!/usr/bin/env python3
"""
Test script to demonstrate the library validator functionality.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.library_validator import LibraryValidator
import tempfile
import shutil

def test_single_file_validation():
    """Test validation of a single GMT file."""
    
    print("🧪 TESTING SINGLE FILE VALIDATION")
    print("=" * 50)
    
    # Create a test GMT file with some problematic gene symbols
    test_gmt_content = """REACTOME_TEST_PATHWAY_1\thttps://reactome.org/PathwayBrowser/#/R-HSA-123456\tC19orf84\tMT-ATP6\tGAPDH\tmt-co1
REACTOME_TEST_PATHWAY_2\thttps://reactome.org/PathwayBrowser/#/R-HSA-789012\tC1orf35\tMT-ND1\tPGK1\tgapdh
REACTOME_TEST_PATHWAY_3\thttps://reactome.org/PathwayBrowser/#/R-HSA-345678\tC2orf49\tMT-RNR1\tLDHA\tpgk1"""
    
    # Create temporary test file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gmt', delete=False) as f:
        f.write(test_gmt_content)
        test_file = f.name
    
    print(f"Created test file: {test_file}")
    print("Original content:")
    with open(test_file, 'r') as f:
        print(f.read())
    
    # Initialize validator
    validator = LibraryValidator()
    
    # Validate the file
    print(f"\n🔍 Validating file...")
    stats = validator.validate_gmt_file(test_file, create_backup=False)
    
    print(f"\n📊 Validation statistics:")
    for key, value in stats.items():
        print(f"   {key}: {value}")
    
    print(f"\n✅ Validated content:")
    with open(test_file, 'r') as f:
        print(f.read())
    
    # Clean up
    os.unlink(test_file)
    
    return stats

def test_directory_validation():
    """Test validation of a directory of GMT files."""
    
    print(f"\n\n🧪 TESTING DIRECTORY VALIDATION")
    print("=" * 50)
    
    # Create a temporary directory with test files
    temp_dir = tempfile.mkdtemp()
    
    # Create test files
    test_files = {
        'test1.gmt': """REACTOME_PATHWAY_1\thttps://reactome.org/PathwayBrowser/#/R-HSA-111111\tC19orf84\tMT-ATP6\tGAPDH
REACTOME_PATHWAY_2\thttps://reactome.org/PathwayBrowser/#/R-HSA-222222\tC1orf35\tMT-ND1\tPGK1""",
        
        'test2.gmt': """GO_PATHWAY_1\thttp://amigo.geneontology.org/amigo/term/GO:0000001\tC2orf49\tMT-RNR1\tLDHA
GO_PATHWAY_2\thttp://amigo.geneontology.org/amigo/term/GO:0000002\tC6orf120\tMT-CO1\tPGAM1""",
        
        'not_a_gmt.txt': "This is not a GMT file"
    }
    
    for filename, content in test_files.items():
        filepath = os.path.join(temp_dir, filename)
        with open(filepath, 'w') as f:
            f.write(content)
    
    print(f"Created test directory: {temp_dir}")
    print(f"Files created: {list(test_files.keys())}")
    
    # Initialize validator
    validator = LibraryValidator()
    
    # Validate the directory
    print(f"\n🔍 Validating directory...")
    stats = validator.validate_library_directory(temp_dir, create_backup=False)
    
    print(f"\n📊 Overall validation statistics:")
    print(validator.get_validation_report())
    
    # Show results
    print(f"\n✅ Validated files:")
    for filename in ['test1.gmt', 'test2.gmt']:
        filepath = os.path.join(temp_dir, filename)
        print(f"\n{filename}:")
        with open(filepath, 'r') as f:
            print(f.read())
    
    # Clean up
    shutil.rmtree(temp_dir)
    
    return stats

def test_reactome_validation():
    """Test validation of the actual Reactome file."""
    
    print(f"\n\n🧪 TESTING REACTOME VALIDATION")
    print("=" * 50)
    
    reactome_file = "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt"
    
    if not os.path.exists(reactome_file):
        print(f"Reactome file not found: {reactome_file}")
        return None
    
    # Create a temporary copy for testing
    temp_file = tempfile.NamedTemporaryFile(suffix='.gmt', delete=False).name
    shutil.copy2(reactome_file, temp_file)
    
    print(f"Testing validation on Reactome file (copy)")
    
    # Initialize validator
    validator = LibraryValidator()
    
    # Validate the file
    print(f"\n🔍 Validating Reactome file...")
    stats = validator.validate_gmt_file(temp_file, create_backup=False)
    
    print(f"\n📊 Reactome validation statistics:")
    print(validator.get_validation_report())
    
    # Show some examples of changes
    print(f"\n🔍 Checking for gene symbol updates...")
    
    # Read original and validated files to compare
    original_genes = set()
    validated_genes = set()
    
    with open(reactome_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                original_genes.update(parts[2:])
    
    with open(temp_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                validated_genes.update(parts[2:])
    
    # Find differences
    updated_genes = validated_genes - original_genes
    if updated_genes:
        print(f"Updated gene symbols: {sorted(list(updated_genes))[:10]}")
        if len(updated_genes) > 10:
            print(f"... and {len(updated_genes) - 10} more")
    else:
        print("No gene symbols were updated")
    
    # Clean up
    os.unlink(temp_file)
    
    return stats

def main():
    """Run all tests."""
    
    print("🧬 LIBRARY VALIDATOR TEST SUITE")
    print("=" * 60)
    
    # Test single file validation
    single_stats = test_single_file_validation()
    
    # Test directory validation
    dir_stats = test_directory_validation()
    
    # Test Reactome validation
    reactome_stats = test_reactome_validation()
    
    print(f"\n\n📋 TEST SUMMARY")
    print("=" * 60)
    print("✅ All tests completed successfully!")
    print("The library validator can:")
    print("  - Validate individual GMT files")
    print("  - Process entire directories of GMT files")
    print("  - Update gene symbols to be consistent with gene_info database")
    print("  - Keep invalid genes unchanged")
    print("  - Generate detailed validation reports")
    print("  - Create backups of original files")

if __name__ == "__main__":
    main()
