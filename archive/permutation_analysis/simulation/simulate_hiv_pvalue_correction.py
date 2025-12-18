#!/usr/bin/env python3
"""
Simulation script: Run iGEA on HIV gene list and show p-value correction impact.

This script demonstrates the post-processing approach to p-value correction:
1. Run iGEA normally (computes raw p-values)
2. Apply correction afterward
3. Compare raw vs corrected p-values

Uses the 5 libraries for which we have permutation data:
- GO BP (Gene Ontology Biological Process)
- GO CC (Gene Ontology Cellular Component)
- GO MF (Gene Ontology Molecular Function)
- KEGG (KEGG Pathways)
- Reactome (Reactome Pathways)
"""

import sys
import logging
from pathlib import Path
from typing import Dict, List
import pandas as pd
import numpy as np

# Add code directory to path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "code"))

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from iter_enrichment import IterativeEnrichment
from pvalue_corrector import load_corrector_from_permutations
from gene_converter import GeneConverter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Paths
DATA_DIR = PROJECT_ROOT / "data"
LIBRARIES_DIR = DATA_DIR / "libraries"
BACKGROUNDS_DIR = DATA_DIR / "backgrounds"
GENE_LISTS_DIR = DATA_DIR / "gene_lists"
PERMUTATION_FILE = PROJECT_ROOT / "permutation_results" / "merged_permutation_results.tsv"

# Libraries to use (matching permutation data)
LIBRARIES = {
    "Reactome": "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
    "KEGG": "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",
    "GO BP": "c5.go.bp.v2025.1.Hs.symbols.gmt",
    "GO MF": "c5.go.mf.v2025.1.Hs.symbols.gmt",
    "GO CC": "c5.go.cc.v2025.1.Hs.symbols.gmt",
}

# iGEA parameters (matching permutation analysis)
PARAMS = {
    "p_threshold": 0.05,
    "min_overlap": 3,
    "min_term_size": 10,
    "max_term_size": 600,
    "max_iterations": 10,
    "p_value_method": "Fisher's Exact Test",
}


def load_hiv_gene_list() -> List[str]:
    """Load HIV gene list and convert Entrez IDs to gene symbols."""
    logger.info("Loading HIV gene list...")
    
    hiv_file = GENE_LISTS_DIR / "HIV.InputGeneList.txt"
    if not hiv_file.exists():
        raise FileNotFoundError(f"HIV gene list not found: {hiv_file}")
    
    # Read Entrez IDs
    with open(hiv_file, 'r') as f:
        entrez_ids = [line.strip() for line in f if line.strip()]
    
    logger.info(f"Loaded {len(entrez_ids)} Entrez IDs from HIV gene list")
    
    # Convert to gene symbols
    converter = GeneConverter()
    gene_symbols = []
    failed = []
    
    for entrez_id in entrez_ids:
        symbol = converter.get_symbol(entrez_id)
        if symbol:
            gene_symbols.append(symbol)
        else:
            failed.append(entrez_id)
    
    if failed:
        logger.warning(f"Failed to convert {len(failed)} Entrez IDs to symbols")
        logger.debug(f"Failed IDs: {failed[:10]}...")  # Show first 10
    
    logger.info(f"Converted to {len(gene_symbols)} gene symbols")
    return gene_symbols


def load_libraries() -> Dict[str, GeneSetLibrary]:
    """Load the 5 libraries for which we have permutation data."""
    logger.info("Loading gene set libraries...")
    libraries = {}
    
    for lib_name, filename in LIBRARIES.items():
        lib_path = LIBRARIES_DIR / filename
        if not lib_path.exists():
            logger.warning(f"Library file not found: {lib_path}, skipping {lib_name}")
            continue
        
        logger.info(f"Loading library: {lib_name}")
        library = GeneSetLibrary(str(lib_path), name=lib_name)
        
        # Filter terms by size
        filtered_terms = [
            t for t in library.library
            if PARAMS["min_term_size"] <= t["size"] <= PARAMS["max_term_size"]
        ]
        library.library = filtered_terms
        library.num_terms = len(filtered_terms)
        library.unique_genes = library.compute_unique_genes()
        library.size = len(library.unique_genes)
        
        libraries[lib_name] = library
        logger.info(f"  {lib_name}: {library.num_terms} terms, {library.size} unique genes")
    
    if not libraries:
        raise ValueError("No libraries could be loaded!")
    
    return libraries


def load_background() -> BackgroundGeneSet:
    """Load the background gene set."""
    bg_path = BACKGROUNDS_DIR / "all_genes.txt"
    if not bg_path.exists():
        raise FileNotFoundError(f"Background file not found: {bg_path}")
    
    logger.info(f"Loading background: {bg_path}")
    bg = BackgroundGeneSet(str(bg_path), name="all_genes", input_format="symbols", skip_validation=True)
    logger.info(f"Loaded background: {bg.size} genes")
    return bg


def run_igea_for_library(
    gene_set: GeneSet,
    library: GeneSetLibrary,
    background: BackgroundGeneSet
) -> IterativeEnrichment:
    """Run iGEA for a single library."""
    logger.info(f"Running iGEA for {library.name}...")
    
    iter_enrich = IterativeEnrichment(
        gene_set=gene_set,
        gene_set_library=library,
        background_gene_set=background,
        min_term_size=PARAMS["min_term_size"],
        max_term_size=PARAMS["max_term_size"],
        p_value_method_name=PARAMS["p_value_method"],
        p_threshold=PARAMS["p_threshold"],
        max_iterations=PARAMS["max_iterations"],
        min_overlap=PARAMS["min_overlap"],
        progress_callback=None,
        use_multiprocessing=True,
    )
    
    logger.info(f"  {library.name}: {len(iter_enrich.results)} iterations")
    return iter_enrich


def apply_correction(
    iter_enrich: IterativeEnrichment,
    corrector,
    gene_list_size: int
) -> pd.DataFrame:
    """Apply p-value correction to iGEA results."""
    # Convert to DataFrame
    df = iter_enrich.to_dataframe()
    
    if df.empty:
        logger.warning("No results to correct")
        return df
    
    # Add gene list size (needed for correction)
    df['gene_list_size'] = gene_list_size
    
    # Apply correction
    df_corrected = corrector.correct_dataframe(df)
    
    return df_corrected


def print_comparison(df_raw: pd.DataFrame, df_corrected: pd.DataFrame, library_name: str):
    """Print comparison of raw vs corrected p-values."""
    print("\n" + "=" * 80)
    print(f"Library: {library_name}")
    print("=" * 80)
    
    if df_raw.empty:
        print("No results for this library")
        return
    
    # Merge for comparison
    comparison = df_raw[['Iteration', 'Term', 'iteration p-value']].copy()
    comparison['raw p-value'] = comparison['iteration p-value']
    comparison['corrected p-value'] = df_corrected['corrected p-value']
    comparison['correction_factor'] = comparison['corrected p-value'] / comparison['raw p-value']
    
    # Display results
    print(f"\nTotal iterations: {len(comparison)}")
    print(f"\n{'Iter':<5} {'Term':<50} {'Raw P':<12} {'Corrected P':<15} {'Factor':<10}")
    print("-" * 100)
    
    for _, row in comparison.iterrows():
        term_short = row['Term'][:47] + "..." if len(row['Term']) > 50 else row['Term']
        print(f"{int(row['Iteration']):<5} {term_short:<50} "
              f"{row['raw p-value']:<12.6f} {row['corrected p-value']:<15.6f} "
              f"{row['correction_factor']:<10.2f}x")
    
    # Summary statistics
    print("\n" + "-" * 100)
    print("Summary Statistics:")
    print(f"  Mean raw p-value: {comparison['raw p-value'].mean():.6f}")
    print(f"  Mean corrected p-value: {comparison['corrected p-value'].mean():.6f}")
    print(f"  Mean correction factor: {comparison['correction_factor'].mean():.2f}x")
    print(f"  Median correction factor: {comparison['correction_factor'].median():.2f}x")
    print(f"  Min correction factor: {comparison['correction_factor'].min():.2f}x")
    print(f"  Max correction factor: {comparison['correction_factor'].max():.2f}x")
    
    # Count significant before/after correction (using 0.05 threshold)
    sig_raw = (comparison['raw p-value'] <= 0.05).sum()
    sig_corrected = (comparison['corrected p-value'] <= 0.05).sum()
    print(f"\n  Significant at 0.05 threshold:")
    print(f"    Raw p-values: {sig_raw}/{len(comparison)} ({100*sig_raw/len(comparison):.1f}%)")
    print(f"    Corrected p-values: {sig_corrected}/{len(comparison)} ({100*sig_corrected/len(comparison):.1f}%)")


def main():
    """Main simulation function."""
    print("\n" + "=" * 80)
    print("HIV Gene List - P-Value Correction Simulation")
    print("=" * 80)
    
    # 1. Load HIV gene list
    print("\n[1/5] Loading HIV gene list...")
    hiv_genes = load_hiv_gene_list()
    print(f"✓ Loaded {len(hiv_genes)} genes")
    
    # 2. Load background
    print("\n[2/5] Loading background gene set...")
    background = load_background()
    print(f"✓ Loaded background: {background.size} genes")
    
    # 3. Create gene set
    print("\n[3/5] Creating gene set...")
    gene_set = GeneSet(
        gene_list=hiv_genes,
        validation_set=background.genes,
        hgcn=False,
        format=False
    )
    print(f"✓ Gene set created: {gene_set.size} genes")
    
    # 4. Load libraries
    print("\n[4/5] Loading gene set libraries...")
    libraries = load_libraries()
    print(f"✓ Loaded {len(libraries)} libraries")
    
    # 5. Load p-value corrector
    print("\n[5/5] Loading p-value corrector...")
    if not PERMUTATION_FILE.exists():
        logger.error(f"Permutation file not found: {PERMUTATION_FILE}")
        print("⚠️  Cannot load corrector - will show raw p-values only")
        corrector = None
    else:
        corrector = load_corrector_from_permutations(str(PERMUTATION_FILE))
        print(f"✓ Loaded corrector with {len([c for c in corrector.stratum_cdfs.values() if c is not None])} stratum CDFs")
    
    # 6. Run iGEA for each library
    print("\n" + "=" * 80)
    print("Running iGEA for each library...")
    print("=" * 80)
    
    all_results = {}
    all_results_corrected = {}
    
    for lib_name, library in libraries.items():
        try:
            # Run iGEA
            iter_enrich = run_igea_for_library(gene_set, library, background)
            all_results[lib_name] = iter_enrich
            
            # Apply correction if corrector is available
            if corrector is not None:
                df_corrected = apply_correction(iter_enrich, corrector, gene_set.size)
                all_results_corrected[lib_name] = df_corrected
            else:
                all_results_corrected[lib_name] = None
            
        except Exception as e:
            logger.error(f"Error processing {lib_name}: {e}", exc_info=True)
            all_results[lib_name] = None
            all_results_corrected[lib_name] = None
    
    # 7. Display results
    print("\n" + "=" * 80)
    print("RESULTS: Raw vs Corrected P-Values")
    print("=" * 80)
    
    for lib_name in libraries.keys():
        if all_results[lib_name] is None:
            continue
        
        df_raw = all_results[lib_name].to_dataframe()
        df_corrected = all_results_corrected[lib_name] if all_results_corrected[lib_name] is not None else None
        
        if df_corrected is not None:
            print_comparison(df_raw, df_corrected, lib_name)
        else:
            print("\n" + "=" * 80)
            print(f"Library: {lib_name}")
            print("=" * 80)
            print("(Corrector not available - showing raw p-values only)")
            if not df_raw.empty:
                print(df_raw[['Iteration', 'Term', 'iteration p-value']].to_string(index=False))
    
    # 8. Save results
    print("\n" + "=" * 80)
    print("Saving results...")
    print("=" * 80)
    
    output_dir = PROJECT_ROOT / "results" / "hiv_simulation"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save combined results
    all_dfs = []
    for lib_name, df_corrected in all_results_corrected.items():
        if df_corrected is not None and not df_corrected.empty:
            all_dfs.append(df_corrected)
    
    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        output_file = output_dir / "hiv_igea_with_correction.tsv"
        combined_df.to_csv(output_file, sep="\t", index=False)
        print(f"✓ Saved combined results to: {output_file}")
        print(f"  Total rows: {len(combined_df)}")
    
    # Save raw results for comparison
    all_dfs_raw = []
    for lib_name, iter_enrich in all_results.items():
        if iter_enrich is not None:
            df_raw = iter_enrich.to_dataframe()
            if not df_raw.empty:
                all_dfs_raw.append(df_raw)
    
    if all_dfs_raw:
        combined_df_raw = pd.concat(all_dfs_raw, ignore_index=True)
        output_file_raw = output_dir / "hiv_igea_raw.tsv"
        combined_df_raw.to_csv(output_file_raw, sep="\t", index=False)
        print(f"✓ Saved raw results to: {output_file_raw}")
    
    print("\n" + "=" * 80)
    print("Simulation complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
