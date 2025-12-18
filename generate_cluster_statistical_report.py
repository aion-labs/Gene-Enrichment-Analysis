#!/usr/bin/env python3
"""
Generate a formatted statistical report for cluster-by-cluster analysis.

This script reads the cluster-by-cluster TSV file and generates a human-readable
statistical report showing benchmark results for each cluster.
"""

import sys
import pandas as pd
from pathlib import Path
from typing import Dict, List


def format_status(z_score: float, percentile: float, is_available: bool = True) -> str:
    """
    Format status indicator based on z-score and percentile.
    
    Args:
        z_score: Z-score value
        percentile: Percentile value
        is_available: Whether the metric is available in null distribution
        
    Returns:
        Status string
    """
    if not is_available or (z_score == 0.0 and percentile == 50.0):
        return "Not available"
    
    if z_score > 2.0 and percentile > 95.0:
        return "✓✓ SIGNIFICANTLY BETTER"
    elif z_score > 1.0 and percentile > 84.0:
        return "✓ BETTER"
    elif z_score > -1.0 and percentile > 16.0:
        return "~ SIMILAR"
    elif z_score < -1.0 and percentile < 16.0:
        return "✗ WORSE"
    else:
        return "~ SIMILAR"


def generate_cluster_report(
    cluster_tsv_file: str,
    output_file: str,
    gene_list_name: str = "Gene List"
) -> None:
    """
    Generate a formatted statistical report from cluster TSV file.
    
    Args:
        cluster_tsv_file: Path to cluster-by-cluster TSV file
        output_file: Path to output report file
        gene_list_name: Name of the gene list being analyzed
    """
    # Read cluster data
    df = pd.read_csv(cluster_tsv_file, sep='\t')
    
    # Generate report
    lines = []
    
    # Header
    lines.append("=" * 100)
    lines.append("CLUSTER-BY-CLUSTER STATISTICAL REPORT")
    lines.append("=" * 100)
    lines.append("")
    lines.append(f"Gene List: {gene_list_name}")
    lines.append(f"Total Clusters: {len(df)}")
    lines.append("")
    lines.append("This report shows each cluster with benchmark statistics comparing")
    lines.append("against random gene lists of similar size.")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("  • Z-score > 2.0 and Percentile > 95%: Significantly better than random (top 5%)")
    lines.append("  • Z-score > 1.0 and Percentile > 84%: Better than random (top 16%)")
    lines.append("  • Z-score ~ 0.0 and Percentile ~ 50%: Similar to random")
    lines.append("  • Z-score < -1.0 and Percentile < 16%: Worse than random")
    lines.append("")
    lines.append("Note: Metrics showing Z=0.0 may not be available in the null distribution.")
    lines.append("      The null distribution may need to be rebuilt with cluster-specific metrics.")
    lines.append("")
    lines.append("=" * 100)
    lines.append("")
    
    # Process each cluster
    for idx, row in df.iterrows():
        cluster_num = int(row['Cluster_Number'])
        
        lines.append("-" * 100)
        lines.append(f"CLUSTER {cluster_num} (Ranked by Size)")
        lines.append("-" * 100)
        lines.append("")
        
        # Basic metrics
        lines.append("Basic Cluster Metrics:")
        lines.append(f"  Cluster Size:        {int(row['Cluster_Size'])} nodes (genes + terms)")
        lines.append(f"  Number of Genes:     {int(row['N_Genes'])}")
        lines.append(f"  Number of Terms:     {int(row['N_Terms'])}")
        lines.append(f"  Number of Edges:     {int(row['N_Edges'])}")
        lines.append(f"  Cluster Density:     {row['Cluster_Density']:.4f}")
        lines.append(f"  Number of Libraries: {int(row['N_Libraries'])}")
        lines.append(f"  Libraries:           {row['Libraries']}")
        lines.append("")
        
        # Statistical benchmarks
        lines.append("Statistical Benchmarks vs Random Gene Lists:")
        lines.append("  Metric                    Value      Z-score   Percentile   Status")
        lines.append("  " + "-" * 90)
        
        # Cluster Size
        size_z = row.get('Cluster_Size_Z', 0.0)
        size_pct = row.get('Cluster_Size_Pct', 50.0)
        size_available = not (size_z == 0.0 and size_pct == 50.0)
        size_status = format_status(size_z, size_pct, size_available)
        lines.append(f"  Cluster Size              {int(row['Cluster_Size']):>8}  {size_z:>8.2f}  {size_pct:>8.1f}%  {size_status}")
        
        # Number of Genes
        genes_z = row.get('N_Genes_Z', 0.0)
        genes_pct = row.get('N_Genes_Pct', 50.0)
        genes_available = not (genes_z == 0.0 and genes_pct == 50.0)
        genes_status = format_status(genes_z, genes_pct, genes_available)
        lines.append(f"  Number of Genes           {int(row['N_Genes']):>8}  {genes_z:>8.2f}  {genes_pct:>8.1f}%  {genes_status}")
        
        # Number of Terms
        terms_z = row.get('N_Terms_Z', 0.0)
        terms_pct = row.get('N_Terms_Pct', 50.0)
        terms_available = not (terms_z == 0.0 and terms_pct == 50.0)
        terms_status = format_status(terms_z, terms_pct, terms_available)
        lines.append(f"  Number of Terms           {int(row['N_Terms']):>8}  {terms_z:>8.2f}  {terms_pct:>8.1f}%  {terms_status}")
        
        # Number of Edges
        edges_z = row.get('N_Edges_Z', 0.0)
        edges_pct = row.get('N_Edges_Pct', 50.0)
        edges_available = not (edges_z == 0.0 and edges_pct == 50.0)
        edges_status = format_status(edges_z, edges_pct, edges_available)
        lines.append(f"  Number of Edges           {int(row['N_Edges']):>8}  {edges_z:>8.2f}  {edges_pct:>8.1f}%  {edges_status}")
        
        # Cluster Density
        density_z = row.get('Cluster_Density_Z', 0.0)
        density_pct = row.get('Cluster_Density_Pct', 50.0)
        density_available = not (density_z == 0.0 and density_pct == 50.0)
        density_status = format_status(density_z, density_pct, density_available)
        lines.append(f"  Cluster Density        {row['Cluster_Density']:>8.4f}  {density_z:>8.2f}  {density_pct:>8.1f}%  {density_status}")
        
        # Number of Libraries
        lib_z = row.get('N_Libraries_Z', 0.0)
        lib_pct = row.get('N_Libraries_Pct', 50.0)
        lib_available = not (lib_z == 0.0 and lib_pct == 50.0)
        lib_status = format_status(lib_z, lib_pct, lib_available)
        lines.append(f"  Number of Libraries       {int(row['N_Libraries']):>8}  {lib_z:>8.2f}  {lib_pct:>8.1f}%  {lib_status}")
        
        lines.append("")
        
        # Terms in cluster
        terms_str = row.get('Terms', '')
        if terms_str and pd.notna(terms_str):
            terms_list = [t.strip() for t in terms_str.split(';') if t.strip()]
            n_terms_total = int(row.get('N_Terms_Total', len(terms_list)))
            
            lines.append(f"Terms in Cluster ({n_terms_total} total):")
            for i, term in enumerate(terms_list[:20], 1):  # Show first 20 terms
                lines.append(f"  {i:2d}. {term}")
            if len(terms_list) > 20:
                lines.append(f"  ... and {len(terms_list) - 20} more terms")
            lines.append("")
    
    # Footer
    lines.append("=" * 100)
    lines.append("END OF REPORT")
    lines.append("=" * 100)
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"✓ Generated statistical report: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: generate_cluster_statistical_report.py <cluster_tsv_file> <output_file> [gene_list_name]")
        sys.exit(1)
    
    cluster_file = sys.argv[1]
    output_file = sys.argv[2]
    gene_list_name = sys.argv[3] if len(sys.argv) > 3 else "Gene List"
    
    generate_cluster_report(cluster_file, output_file, gene_list_name)
