#!/usr/bin/env python3
from pathlib import Path
import sys

print("Current working directory:", Path.cwd())
print("Python path:", sys.path[:3])

# Test the path calculation that gene_converter uses
code_file = Path("code/gene_converter.py")
print("Code file exists:", code_file.exists())
print("Code file resolved:", code_file.resolve())
print("Parent:", code_file.resolve().parent)
print("Parent.parent:", code_file.resolve().parent.parent)

# Test the gene info path
gene_info_path = code_file.resolve().parent.parent / "data" / "recent_release" / "Homo_sapiens.gene_info"
print("Gene info path:", gene_info_path)
print("Gene info exists:", gene_info_path.exists())

# Test the correct path
correct_path = Path.cwd() / "data" / "recent_release" / "Homo_sapiens.gene_info"
print("Correct path:", correct_path)
print("Correct path exists:", correct_path.exists())



