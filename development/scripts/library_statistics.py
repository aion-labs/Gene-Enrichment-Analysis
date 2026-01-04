#!/usr/bin/env python3
"""
Comprehensive library statistics generator for gene set libraries.
Analyzes all available libraries and provides detailed statistics.
"""

import sys
import os
import statistics
import json
import re
from pathlib import Path
from collections import Counter
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'code'))

from code.gene_set_library import GeneSetLibrary

def load_aliases():
    """Load library aliases from alias.json file."""
    alias_file = Path("data/libraries/alias.json")
    if alias_file.exists():
        with open(alias_file, 'r') as f:
            return json.load(f)
    return []

def extract_category(filename: str) -> str:
    """Extract category (C1-C8, H) from filename."""
    # Extract the category from the filename
    match = re.match(r'^([ch]\d*|[h])\.', filename.lower())
    if match:
        category = match.group(1).upper()
        # Handle special cases
        if category == 'H':
            return 'H'
        elif category.startswith('C'):
            return category
    return 'Other'

def get_friendly_name(filename: str, aliases: list) -> str:
    """Get friendly name from aliases or generate one from filename."""
    filename_only = Path(filename).name
    
    # Look for exact match in aliases
    for alias in aliases:
        if alias['file'] == filename_only:
            return alias['name']
    
    # If no alias found, generate a friendly name from filename
    category = extract_category(filename_only)
    if category == 'H':
        return 'H: Hallmark Gene Sets'
    elif category.startswith('C'):
        # Extract subcategory if available
        parts = filename_only.split('.')
        if len(parts) >= 2:
            subcategory = parts[1].upper()
            if subcategory in ['ALL', 'GO', 'MIR', 'TFT', 'CGP', 'CP', '3CA', 'CGN', 'CM', 'VAX']:
                return f"{category}: {subcategory}"
        return f"{category}: Gene Sets"
    
    return Path(filename).stem

def analyze_library(library_path: str, aliases: list) -> dict:
    """
    Analyze a single library and return comprehensive statistics.
    
    Args:
        library_path: Path to the .gmt file
        aliases: List of library aliases
        
    Returns:
        Dictionary containing all statistics
    """
    try:
        lib = GeneSetLibrary(library_path, name=Path(library_path).stem)
    except Exception as e:
        print(f"Error loading library {library_path}: {e}")
        return None
    
    filename = Path(library_path).name
    category = extract_category(filename)
    friendly_name = get_friendly_name(library_path, aliases)
    
    # Basic statistics
    stats = {
        'name': lib.name,
        'friendly_name': friendly_name,
        'category': category,
        'file_path': library_path,
        'total_terms': lib.num_terms,
        'unique_genes': lib.size,
        'organism': lib.organism
    }
    
    # Term size statistics
    term_sizes = [term["size"] for term in lib.library]
    stats['term_sizes'] = {
        'min': min(term_sizes),
        'max': max(term_sizes),
        'mean': statistics.mean(term_sizes),
        'median': statistics.median(term_sizes),
        'std': statistics.stdev(term_sizes) if len(term_sizes) > 1 else 0
    }
    
    # Gene frequency analysis
    gene_term_counts = Counter()
    for term in lib.library:
        for gene in term["genes"]:
            gene_term_counts[gene] += 1
    
    stats['gene_frequency'] = {
        'appearing_in_1_term': sum(1 for count in gene_term_counts.values() if count == 1),
        'appearing_in_2_5_terms': sum(1 for count in gene_term_counts.values() if 2 <= count <= 5),
        'appearing_in_6_10_terms': sum(1 for count in gene_term_counts.values() if 6 <= count <= 10),
        'appearing_in_11_20_terms': sum(1 for count in gene_term_counts.values() if 11 <= count <= 20),
        'appearing_in_21_50_terms': sum(1 for count in gene_term_counts.values() if 21 <= count <= 50),
        'appearing_in_50_plus_terms': sum(1 for count in gene_term_counts.values() if count > 50)
    }
    
    # Most common genes
    stats['top_genes'] = gene_term_counts.most_common(10)
    
    # Term size distribution
    size_ranges = {
        'small (1-10)': sum(1 for size in term_sizes if 1 <= size <= 10),
        'medium (11-50)': sum(1 for size in term_sizes if 11 <= size <= 50),
        'large (51-200)': sum(1 for size in term_sizes if 51 <= size <= 200),
        'very_large (201-500)': sum(1 for size in term_sizes if 201 <= size <= 500),
        'huge (500+)': sum(1 for size in term_sizes if size > 500)
    }
    stats['term_size_distribution'] = size_ranges
    
    return stats

def print_library_stats(stats: dict, detailed: bool = False):
    """Print formatted statistics for a library."""
    if not stats:
        return
    
    print(f"\n📚 {stats['friendly_name']}")
    print(f"🏷️  Category: {stats['category']}")
    print("=" * 60)
    print(f"📊 Total terms: {stats['total_terms']:,}")
    print(f"🧬 Unique genes: {stats['unique_genes']:,}")
    print(f"🦠 Organism: {stats['organism']}")
    
    # Term size statistics
    sizes = stats['term_sizes']
    print(f"\n📏 Term Size Statistics:")
    print(f"   Smallest term: {sizes['min']} genes")
    print(f"   Largest term: {sizes['max']} genes")
    print(f"   Average term size: {sizes['mean']:.1f} genes")
    print(f"   Median term size: {sizes['median']:.1f} genes")
    print(f"   Standard deviation: {sizes['std']:.1f} genes")
    
    # Term size distribution
    dist = stats['term_size_distribution']
    print(f"\n📊 Term Size Distribution:")
    for range_name, count in dist.items():
        percentage = (count / stats['total_terms']) * 100
        print(f"   {range_name}: {count:,} terms ({percentage:.1f}%)")
    
    # Gene frequency
    freq = stats['gene_frequency']
    print(f"\n📈 Gene Frequency Analysis:")
    print(f"   Genes appearing in 1 term: {freq['appearing_in_1_term']:,}")
    print(f"   Genes appearing in 2-5 terms: {freq['appearing_in_2_5_terms']:,}")
    print(f"   Genes appearing in 6-10 terms: {freq['appearing_in_6_10_terms']:,}")
    print(f"   Genes appearing in 11-20 terms: {freq['appearing_in_11_20_terms']:,}")
    print(f"   Genes appearing in 21-50 terms: {freq['appearing_in_21_50_terms']:,}")
    print(f"   Genes appearing in 50+ terms: {freq['appearing_in_50_plus_terms']:,}")
    
    if detailed and stats['top_genes']:
        print(f"\n🏆 Top 10 most common genes:")
        for gene, count in stats['top_genes']:
            print(f"   {gene}: {count} terms")

def generate_summary_table(all_stats: list) -> pd.DataFrame:
    """Generate a summary table of all libraries."""
    summary_data = []
    
    for stats in all_stats:
        if not stats:
            continue
            
        summary_data.append({
            'Category': stats['category'],
            'Library': stats['friendly_name'],
            'Original Name': stats['name'],
            'Total Terms': stats['total_terms'],
            'Unique Genes': stats['unique_genes'],
            'Avg Term Size': round(stats['term_sizes']['mean'], 1),
            'Median Term Size': round(stats['term_sizes']['median'], 1),
            'Min Term Size': stats['term_sizes']['min'],
            'Max Term Size': stats['term_sizes']['max'],
            'Genes in 1 Term': stats['gene_frequency']['appearing_in_1_term'],
            'Genes in 50+ Terms': stats['gene_frequency']['appearing_in_50_plus_terms']
        })
    
    return pd.DataFrame(summary_data)

def main():
    """Main function to analyze all libraries."""
    print("🔬 COMPREHENSIVE LIBRARY STATISTICS")
    print("=" * 80)
    
    # Load aliases
    aliases = load_aliases()
    print(f"📋 Loaded {len(aliases)} library aliases")
    
    # Get all .gmt files from the libraries directory
    libraries_dir = Path("data/libraries")
    gmt_files = list(libraries_dir.glob("*.gmt"))
    
    if not gmt_files:
        print("❌ No .gmt files found in data/libraries/")
        return
    
    print(f"📁 Found {len(gmt_files)} library files")
    
    # Analyze each library
    all_stats = []
    for gmt_file in sorted(gmt_files):
        print(f"\n🔄 Analyzing {gmt_file.name}...")
        stats = analyze_library(str(gmt_file), aliases)
        all_stats.append(stats)
        
        if stats:
            print_library_stats(stats, detailed=False)
    
    # Generate summary table
    print(f"\n📋 SUMMARY TABLE")
    print("=" * 80)
    summary_df = generate_summary_table(all_stats)
    print(summary_df.to_string(index=False))
    
    # Save detailed results to CSV
    output_file = "library_statistics_summary.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"\n💾 Summary saved to: {output_file}")
    
    # Generate detailed report
    print(f"\n📄 DETAILED REPORT")
    print("=" * 80)
    for stats in all_stats:
        if stats:
            print_library_stats(stats, detailed=True)
    
    # Overall statistics
    print(f"\n🎯 OVERALL STATISTICS")
    print("=" * 80)
    total_terms = sum(stats['total_terms'] for stats in all_stats if stats)
    total_genes = sum(stats['unique_genes'] for stats in all_stats if stats)
    avg_terms_per_lib = total_terms / len([s for s in all_stats if s])
    avg_genes_per_lib = total_genes / len([s for s in all_stats if s])
    
    print(f"📊 Total libraries analyzed: {len([s for s in all_stats if s])}")
    print(f"📊 Total terms across all libraries: {total_terms:,}")
    print(f"🧬 Total unique genes across all libraries: {total_genes:,}")
    print(f"📊 Average terms per library: {avg_terms_per_lib:.1f}")
    print(f"🧬 Average genes per library: {avg_genes_per_lib:.1f}")

if __name__ == "__main__":
    main()
