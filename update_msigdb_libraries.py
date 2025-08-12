#!/usr/bin/env python3
"""
Script to update MSigDB libraries for gene symbol consistency.
This script demonstrates how to use the LibraryValidator to ensure
consistency between gene_info database and MSigDB libraries.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.library_validator import LibraryValidator
import argparse
from pathlib import Path

def update_libraries(library_dir: str, output_dir: str = None, create_backup: bool = True):
    """
    Update MSigDB libraries for gene symbol consistency.
    
    Args:
        library_dir: Directory containing MSigDB GMT files
        output_dir: Output directory (if None, overwrites original files)
        create_backup: Whether to create backups of original files
    """
    
    print("🧬 MSIGDB LIBRARY UPDATE")
    print("=" * 60)
    print(f"Library directory: {library_dir}")
    print(f"Output directory: {output_dir if output_dir else '(overwrite original)'}")
    print(f"Create backups: {create_backup}")
    print()
    
    # Initialize validator
    print("🔧 Initializing library validator...")
    validator = LibraryValidator()
    
    # Validate all libraries
    print("🔍 Validating and updating libraries...")
    stats = validator.validate_library_directory(
        library_dir, 
        output_dir, 
        file_pattern="*.gmt",
        create_backup=create_backup
    )
    
    # Print detailed report
    print("\n" + validator.get_validation_report())
    
    # Additional analysis
    print("\n📊 DETAILED ANALYSIS")
    print("-" * 40)
    
    total_genes = stats['total_genes']
    validated_genes = stats['validated_genes']
    unchanged_genes = stats['unchanged_genes']
    invalid_genes = stats['invalid_genes']
    
    print(f"Total genes processed: {total_genes:,}")
    print(f"Genes updated: {validated_genes:,} ({validated_genes/total_genes*100:.2f}%)")
    print(f"Genes unchanged: {unchanged_genes:,} ({unchanged_genes/total_genes*100:.2f}%)")
    print(f"Genes that could not be validated: {invalid_genes:,} ({invalid_genes/total_genes*100:.2f}%)")
    
    if invalid_genes > 0:
        print(f"\n⚠️  Note: {invalid_genes:,} genes could not be validated.")
        print("   These genes will remain unchanged in the libraries.")
        print("   They will be ignored during enrichment analysis since the")
        print("   background is library-specific and will exclude them anyway.")
    
    print(f"\n✅ Library update completed successfully!")
    
    if output_dir:
        print(f"Updated libraries saved to: {output_dir}")
    else:
        print("Original libraries have been updated in place.")
        if create_backup:
            print("Backup files have been created with timestamp suffixes.")

def main():
    """Main function for command-line interface."""
    
    parser = argparse.ArgumentParser(
        description="Update MSigDB libraries for gene symbol consistency",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Update libraries in place with backups
  python update_msigdb_libraries.py data/libraries/
  
  # Update libraries and save to new directory
  python update_msigdb_libraries.py data/libraries/ -o data/libraries_validated/
  
  # Update without creating backups
  python update_msigdb_libraries.py data/libraries/ --no-backup
        """
    )
    
    parser.add_argument(
        "library_dir",
        help="Directory containing MSigDB GMT files"
    )
    
    parser.add_argument(
        "-o", "--output",
        help="Output directory for updated libraries (if not specified, overwrites original)"
    )
    
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Don't create backup files when overwriting original files"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose logging"
    )
    
    args = parser.parse_args()
    
    # Validate input directory
    library_path = Path(args.library_dir)
    if not library_path.exists():
        print(f"❌ Error: Library directory not found: {args.library_dir}")
        return 1
    
    if not library_path.is_dir():
        print(f"❌ Error: Path is not a directory: {args.library_dir}")
        return 1
    
    # Check for GMT files
    gmt_files = list(library_path.glob("*.gmt"))
    if not gmt_files:
        print(f"❌ Error: No GMT files found in {args.library_dir}")
        return 1
    
    print(f"Found {len(gmt_files)} GMT files to process:")
    for gmt_file in gmt_files:
        print(f"  - {gmt_file.name}")
    print()
    
    # Set up logging
    import logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')
    
    try:
        # Update libraries
        update_libraries(
            args.library_dir,
            args.output,
            create_backup=not args.no_backup
        )
        return 0
        
    except Exception as e:
        print(f"❌ Error updating libraries: {e}")
        return 1

if __name__ == "__main__":
    exit(main())
