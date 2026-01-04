#!/usr/bin/env python3
"""
Example: Network Connectivity Benchmarking

This script demonstrates how to:
1. Compute connectivity metrics from permutation results
2. Build null distributions
3. Benchmark real iGEA results against the null
"""

import sys
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from network_connectivity_benchmark import (
    compute_connectivity_from_permutation_results,
    build_null_distribution,
    benchmark_real_results,
    NetworkConnectivityAnalyzer
)
from iter_enrichment import IterativeEnrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from background_gene_set import BackgroundGeneSet


def example_compute_null_distribution():
    """Example: Compute null distribution from permutation results."""
    print("=" * 60)
    print("Example 1: Computing Null Distribution from Permutations")
    print("=" * 60)
    
    permutation_file = (
        PROJECT_ROOT / "archive" / "permutation_analysis" / "results" /
        "permutation_results" / "merged_permutation_results.tsv"
    )
    
    if not permutation_file.exists():
        print(f"⚠️  Permutation file not found: {permutation_file}")
        print("   Skipping null distribution computation...")
        return None
    
    output_metrics = PROJECT_ROOT / "results" / "connectivity_metrics.tsv"
    output_null = PROJECT_ROOT / "results" / "connectivity_null_distribution.json"
    
    print(f"\n1. Computing connectivity metrics from permutations...")
    print(f"   Input: {permutation_file}")
    
    metrics_df = compute_connectivity_from_permutation_results(
        str(permutation_file),
        str(output_metrics)
    )
    
    print(f"   ✓ Processed {len(metrics_df)} permutations")
    print(f"   ✓ Saved metrics to: {output_metrics}")
    
    print(f"\n2. Building null distribution...")
    null_dist = build_null_distribution(
        metrics_df,
        str(output_null)
    )
    
    print(f"   ✓ Built null distribution for {len(null_dist)} gene list sizes")
    print(f"   ✓ Saved to: {output_null}")
    
    # Show sample statistics
    if null_dist:
        sample_size = list(null_dist.keys())[0]
        sample_stats = null_dist[sample_size]
        print(f"\n3. Sample statistics (gene list size {sample_size}):")
        print(f"   Number of permutations: {sample_stats['n_permutations']}")
        print(f"   Avg connections per gene:")
        print(f"     Mean: {sample_stats['avg_connections_per_gene']['mean']:.3f}")
        print(f"     Std:  {sample_stats['avg_connections_per_gene']['std']:.3f}")
        print(f"     Q5-Q95: {sample_stats['avg_connections_per_gene']['q5']:.3f} - {sample_stats['avg_connections_per_gene']['q95']:.3f}")
    
    return null_dist


def example_benchmark_real_results(null_distribution):
    """Example: Benchmark real iGEA results."""
    print("\n" + "=" * 60)
    print("Example 2: Benchmarking Real iGEA Results")
    print("=" * 60)
    
    if null_distribution is None:
        print("⚠️  No null distribution available. Skipping benchmark...")
        return
    
    # Example: Load a real iGEA result
    # (In practice, you would run iGEA first)
    print("\n1. Loading real iGEA results...")
    print("   (This is a placeholder - replace with your actual iGEA results)")
    
    # Example iGEA results structure
    example_results = [
        {
            'Term': 'GOBP: CELL_CYCLE',
            'Iteration': 1,
            'Genes removed for next iteration': ['CDK1', 'CDK2', 'CCNB1', 'CCNA2']
        },
        {
            'Term': 'GOBP: DNA_REPLICATION',
            'Iteration': 2,
            'Genes removed for next iteration': ['PCNA', 'MCM2', 'MCM3']
        },
    ]
    
    gene_list_size = 50  # Example size
    
    print(f"   Gene list size: {gene_list_size}")
    print(f"   Number of iterations: {len(example_results)}")
    
    print("\n2. Computing connectivity metrics...")
    analyzer = NetworkConnectivityAnalyzer()
    analyzer.add_igea_results(example_results)
    real_metrics = analyzer.compute_metrics()
    
    print(f"   Genes: {real_metrics['n_genes']}")
    print(f"   Terms: {real_metrics['n_terms']}")
    print(f"   Edges: {real_metrics['n_edges']}")
    print(f"   Avg connections per gene: {real_metrics['avg_connections_per_gene']:.2f}")
    
    print("\n3. Benchmarking against null distribution...")
    benchmark = benchmark_real_results(
        example_results,
        gene_list_size,
        null_distribution
    )
    
    if benchmark and 'comparison' in benchmark:
        print("\n4. Benchmark Results:")
        print(f"   {'Metric':<30} {'Real':<10} {'Null Mean':<12} {'Z-score':<10} {'Percentile':<12}")
        print("   " + "-" * 74)
        
        for metric, comp in benchmark['comparison'].items():
            print(f"   {metric:<30} "
                  f"{comp['real_value']:<10.3f} "
                  f"{comp['null_mean']:<12.3f} "
                  f"{comp['z_score']:<10.2f} "
                  f"{comp['percentile']:<12.1f}%")
        
        print("\n5. Interpretation:")
        better_metrics = sum(1 for comp in benchmark['comparison'].values() if comp['is_better'])
        total_metrics = len(benchmark['comparison'])
        print(f"   {better_metrics}/{total_metrics} metrics are better than null")
        
        if better_metrics == total_metrics:
            print("   ✓ Network connectivity is significantly better than random!")
        elif better_metrics > total_metrics / 2:
            print("   ⚠ Network connectivity is moderately better than random")
        else:
            print("   ✗ Network connectivity is not significantly better than random")


def main():
    """Run examples."""
    print("\n" + "=" * 60)
    print("Network Connectivity Benchmarking Examples")
    print("=" * 60)
    
    # Example 1: Compute null distribution
    null_dist = example_compute_null_distribution()
    
    # Example 2: Benchmark real results
    example_benchmark_real_results(null_dist)
    
    print("\n" + "=" * 60)
    print("Examples complete!")
    print("=" * 60)
    print("\nNext steps:")
    print("1. Run the full benchmark on all permutation data")
    print("2. Test on real gene lists (e.g., HIV gene list)")
    print("3. Integrate into iGEA workflow")


if __name__ == "__main__":
    main()
