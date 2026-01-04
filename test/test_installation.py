#!/usr/bin/env python3
"""
Comprehensive installation test script.

This script verifies that the installation is complete and functional by:
1. Checking Python version
2. Checking all required dependencies
3. Testing imports
4. Verifying file structure
5. Running basic CLI and Streamlit tests
"""

import sys
import subprocess
from pathlib import Path

def check_python_version():
    """Check that Python version meets requirements."""
    print("\n" + "=" * 60)
    print("Checking Python Version")
    print("=" * 60)
    
    version = sys.version_info
    print(f"Python version: {version.major}.{version.minor}.{version.micro}")
    
    if version.major == 3 and version.minor >= 12:
        print("✅ Python version meets requirement (>= 3.12)")
        return True
    else:
        print(f"❌ Python version {version.major}.{version.minor} does not meet requirement (>= 3.12)")
        return False

def check_dependencies():
    """Check that all required dependencies are installed."""
    print("\n" + "=" * 60)
    print("Checking Dependencies")
    print("=" * 60)
    
    required_packages = [
        "pandas",
        "scipy",
        "statsmodels",
        "plotly",
        "pydot",
        "networkx",
        "streamlit",
        "typer",
        "PIL",
    ]
    
    all_installed = True
    for package in required_packages:
        try:
            if package == "PIL":
                __import__("PIL")
                print(f"✅ {package} (Pillow) installed")
            else:
                mod = __import__(package)
                version = getattr(mod, "__version__", "unknown")
                print(f"✅ {package} installed (version: {version})")
        except ImportError:
            print(f"❌ {package} not installed")
            all_installed = False
    
    return all_installed

def check_file_structure():
    """Check that the project structure is correct."""
    print("\n" + "=" * 60)
    print("Checking File Structure")
    print("=" * 60)
    
    project_root = Path(__file__).resolve().parent.parent
    required_paths = [
        "code/cli.py",
        "code/streamlit_app.py",
        "data/libraries",
        "data/backgrounds",
        "data/gene_lists",
    ]
    
    all_exist = True
    for path_str in required_paths:
        path = project_root / path_str
        if path.exists():
            print(f"✅ {path_str} exists")
        else:
            print(f"❌ {path_str} not found")
            all_exist = False
    
    return all_exist

def run_cli_test():
    """Run the CLI test script."""
    print("\n" + "=" * 60)
    print("Running CLI Test")
    print("=" * 60)
    
    test_script = Path(__file__).parent / "test_cli.py"
    if not test_script.exists():
        print(f"⚠️  CLI test script not found: {test_script}")
        return False
    
    try:
        result = subprocess.run(
            [sys.executable, str(test_script)],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print("❌ CLI test timed out")
        return False
    except Exception as e:
        print(f"❌ Error running CLI test: {e}")
        return False

def run_streamlit_test():
    """Run the Streamlit test script."""
    print("\n" + "=" * 60)
    print("Running Streamlit Test")
    print("=" * 60)
    
    test_script = Path(__file__).parent / "test_streamlit.py"
    if not test_script.exists():
        print(f"⚠️  Streamlit test script not found: {test_script}")
        return False
    
    try:
        result = subprocess.run(
            [sys.executable, str(test_script)],
            capture_output=True,
            text=True
        )
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        return result.returncode == 0
    except Exception as e:
        print(f"❌ Error running Streamlit test: {e}")
        return False

def main():
    """Run all installation tests."""
    print("\n" + "🔍 " * 20)
    print("Installation Verification Test Suite")
    print("🔍 " * 20 + "\n")
    
    results = {}
    
    # Test 1: Python version
    results["python_version"] = check_python_version()
    
    # Test 2: Dependencies
    results["dependencies"] = check_dependencies()
    
    # Test 3: File structure
    results["file_structure"] = check_file_structure()
    
    # Test 4: CLI
    if results["file_structure"]:
        results["cli"] = run_cli_test()
    else:
        print("\n⚠️  Skipping CLI test (file structure check failed)")
        results["cli"] = False
    
    # Test 5: Streamlit
    if results["dependencies"]:
        results["streamlit"] = run_streamlit_test()
    else:
        print("\n⚠️  Skipping Streamlit test (dependencies check failed)")
        results["streamlit"] = False
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{test_name:20s}: {status}")
    
    all_passed = all(results.values())
    print("=" * 60)
    if all_passed:
        print("✅ ALL INSTALLATION TESTS PASSED!")
        print("   The installation is complete and functional.")
    else:
        print("❌ SOME INSTALLATION TESTS FAILED")
        print("   Please review the output above and fix any issues.")
    print("=" * 60 + "\n")
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())

