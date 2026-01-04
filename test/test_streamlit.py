#!/usr/bin/env python3
"""
Test script to verify Streamlit app installation and basic functionality.

This script checks that:
1. Streamlit can be imported
2. The streamlit_app module can be imported
3. Core dependencies are available
4. Data files are accessible
"""

import sys
from pathlib import Path

def test_imports():
    """Test that all required modules can be imported."""
    print("\n" + "=" * 60)
    print("Testing Imports")
    print("=" * 60)
    
    try:
        import streamlit
        print(f"✅ streamlit imported (version: {streamlit.__version__})")
    except ImportError as e:
        print(f"❌ Failed to import streamlit: {e}")
        return False
    
    try:
        import pandas
        print(f"✅ pandas imported (version: {pandas.__version__})")
    except ImportError as e:
        print(f"❌ Failed to import pandas: {e}")
        return False
    
    try:
        import plotly
        print(f"✅ plotly imported")
    except ImportError as e:
        print(f"❌ Failed to import plotly: {e}")
        return False
    
    # Add code directory to path
    project_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(project_root / "code"))
    
    try:
        from streamlit_app import ROOT
        print(f"✅ streamlit_app imported successfully")
        print(f"   ROOT path: {ROOT}")
        return True
    except ImportError as e:
        print(f"❌ Failed to import streamlit_app: {e}")
        import traceback
        traceback.print_exc()
        return False
    except Exception as e:
        print(f"⚠️  streamlit_app imported but error occurred: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_data_files():
    """Test that required data files exist."""
    print("\n" + "=" * 60)
    print("Testing Data Files")
    print("=" * 60)
    
    project_root = Path(__file__).resolve().parent.parent
    
    # Check data directories
    data_dir = project_root / "data"
    if not data_dir.exists():
        print(f"❌ Data directory not found: {data_dir}")
        return False
    print(f"✅ Data directory exists: {data_dir}")
    
    # Check backgrounds
    backgrounds_dir = data_dir / "backgrounds"
    if backgrounds_dir.exists():
        bg_files = list(backgrounds_dir.glob("*.txt"))
        print(f"✅ Backgrounds directory exists with {len(bg_files)} files")
    else:
        print(f"⚠️  Backgrounds directory not found: {backgrounds_dir}")
    
    # Check libraries
    libraries_dir = data_dir / "libraries"
    if libraries_dir.exists():
        gmt_files = list(libraries_dir.glob("*.gmt"))
        print(f"✅ Libraries directory exists with {len(gmt_files)} GMT files")
    else:
        print(f"❌ Libraries directory not found: {libraries_dir}")
        return False
    
    # Check gene lists
    gene_lists_dir = data_dir / "gene_lists"
    if gene_lists_dir.exists():
        gene_list_files = list(gene_lists_dir.glob("*.txt"))
        print(f"✅ Gene lists directory exists with {len(gene_list_files)} files")
    else:
        print(f"⚠️  Gene lists directory not found: {gene_lists_dir}")
    
    # Check alias.json
    alias_file = libraries_dir / "alias.json"
    if alias_file.exists():
        print(f"✅ Library alias.json exists")
    else:
        print(f"⚠️  Library alias.json not found")
    
    return True

def test_static_files():
    """Test that static files (logos) exist."""
    print("\n" + "=" * 60)
    print("Testing Static Files")
    print("=" * 60)
    
    project_root = Path(__file__).resolve().parent.parent
    static_dir = project_root / "code" / "static"
    
    if static_dir.exists():
        logo_files = list(static_dir.glob("*.png"))
        print(f"✅ Static directory exists with {len(logo_files)} image files")
        return True
    else:
        print(f"⚠️  Static directory not found: {static_dir}")
        return False

def main():
    """Run all Streamlit installation tests."""
    print("\n" + "🧪 " * 20)
    print("Streamlit Installation Test Suite")
    print("🧪 " * 20 + "\n")
    
    all_passed = True
    
    # Test 1: Imports
    if not test_imports():
        all_passed = False
    
    # Test 2: Data files
    if not test_data_files():
        all_passed = False
    
    # Test 3: Static files
    if not test_static_files():
        all_passed = False
    
    print("\n" + "=" * 60)
    if all_passed:
        print("✅ ALL STREAMLIT TESTS PASSED!")
        print("   The Streamlit app should be ready to run.")
        print("   Try: streamlit run code/streamlit_app.py")
    else:
        print("⚠️  SOME TESTS FAILED OR HAD WARNINGS")
        print("   Please check the output above for details.")
    print("=" * 60 + "\n")
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())

