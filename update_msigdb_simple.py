#!/usr/bin/env python3
"""
Simple MSigDB Update Script

This script updates the MSigDB gene set libraries using the zip file
in the data/recent_release directory.

Usage:
    python update_msigdb_simple.py
"""

import sys
from pathlib import Path

# Add the current directory to the path so we can import the main update script
sys.path.append(str(Path(__file__).parent))

from update_msigdb import main as update_main

def main():
    """Simple wrapper that uses the zip file in recent_release."""
    
    # Find the zip file
    zip_file = Path("data/recent_release/msigdb_v2025.1.Hs_files_to_download_locally.zip")
    
    if not zip_file.exists():
        print(f"Error: Zip file not found at {zip_file}")
        print("Please ensure the MSigDB zip file is in the data/recent_release directory")
        return 1
    
    print(f"Found MSigDB zip file: {zip_file}")
    print("Starting update process...")
    print("This will update C2 and C5 gene set libraries, download NCBI gene information files,")
    print("and create a comprehensive background gene list with validated genes.")
    
    # Set up the arguments for the main update script
    sys.argv = [
        "update_msigdb.py",
        "--zip-file",
        str(zip_file)
    ]
    
    # Run the main update script
    return update_main()

if __name__ == "__main__":
    exit(main())
