#!/usr/bin/env python3
"""
Generate permutation distribution of p-values for iGEA analysis.

This script generates random gene lists and runs iterative enrichment analysis (iGEA)
to create a null distribution of p-values for statistical validation.

Usage:
    python code/generate_permutation_distribution.py [--n-permutations N] [--n-jobs N] [--resume]
    
    Or from the code/ directory:
        cd code
        python generate_permutation_distribution.py [options]
    
    Example for c6i.8xlarge (32 vCPUs):
        python code/generate_permutation_distribution.py --n-jobs 32
    
    Results are saved to: results/permutation_results/
"""

import argparse
import json
import logging
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from multiprocessing import Pool, cpu_count, set_start_method
from functools import partial

# Set multiprocessing start method to 'spawn' to allow nested multiprocessing
# This is needed because enrichment code uses multiprocessing internally
try:
    set_start_method('spawn', force=True)
except RuntimeError:
    # Already set, ignore
    pass

# Configuration - define paths first
# When script is in code/ folder, go up one level to project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
LIBRARIES_DIR = DATA_DIR / "libraries"
BACKGROUNDS_DIR = DATA_DIR / "backgrounds"
RESULTS_DIR = PROJECT_ROOT / "results"
OUTPUT_DIR = RESULTS_DIR / "permutation_results"

# Add current directory (code/) to path for imports
# When script is in code/ folder, imports work directly
sys.path.insert(0, str(Path(__file__).parent))

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from iter_enrichment import IterativeEnrichment

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(str(PROJECT_ROOT / "permutation_distribution.log")),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Library mapping (name -> filename)
LIBRARIES = {
    "Reactome": "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
    "KEGG": "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",
    "GO BP": "c5.go.bp.v2025.1.Hs.symbols.gmt",
    "GO MF": "c5.go.mf.v2025.1.Hs.symbols.gmt",
    "GO CC": "c5.go.cc.v2025.1.Hs.symbols.gmt",
}

# Default parameters (matching Streamlit defaults)
DEFAULT_PARAMS = {
    "p_threshold": 0.05,  # User specified
    "min_overlap": 3,
    "min_term_size": 10,
    "max_term_size": 600,
    "max_iterations": 30,  # Increased for permutation analysis
    "p_value_method": "Fisher's Exact Test",
}

# Gene list sizes: 50 to 1000 in increments of 50
GENE_LIST_SIZES = list(range(50, 1001, 50))  # [50, 100, 150, ..., 1000]


def load_background() -> BackgroundGeneSet:
    """Load the all_genes background gene set."""
    bg_path = BACKGROUNDS_DIR / "all_genes.txt"
    if not bg_path.exists():
        raise FileNotFoundError(f"Background file not found: {bg_path}")
    
    logger.info(f"Loading background: {bg_path}")
    bg = BackgroundGeneSet(str(bg_path), name="all_genes", input_format="symbols")
    logger.info(f"Loaded background: {bg.size} genes")
    return bg


def load_libraries() -> Dict[str, GeneSetLibrary]:
    """Load all required gene set libraries."""
    libraries = {}
    
    for lib_name, filename in LIBRARIES.items():
        lib_path = LIBRARIES_DIR / filename
        if not lib_path.exists():
            logger.warning(f"Library file not found: {lib_path}, skipping {lib_name}")
            continue
        
        logger.info(f"Loading library: {lib_name} ({filename})")
        library = GeneSetLibrary(str(lib_path), name=lib_name)
        
        # Filter terms by size (same as Streamlit app)
        filtered_terms = [
            t for t in library.library
            if DEFAULT_PARAMS["min_term_size"] <= t["size"] <= DEFAULT_PARAMS["max_term_size"]
        ]
        library.library = filtered_terms
        library.num_terms = len(filtered_terms)
        library.unique_genes = library.compute_unique_genes()
        library.size = len(library.unique_genes)
        
        libraries[lib_name] = library
        logger.info(f"  Loaded {lib_name}: {library.num_terms} terms, {library.size} unique genes")
    
    if not libraries:
        raise ValueError("No libraries could be loaded!")
    
    return libraries


def create_random_gene_set(
    background: BackgroundGeneSet,
    size: int,
    seed: Optional[int] = None
) -> GeneSet:
    """
    Create a random gene set by sampling from the background.
    
    Args:
        background: Background gene set to sample from
        size: Number of genes to sample
        seed: Random seed for reproducibility
        
    Returns:
        GeneSet object with randomly sampled genes
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Sample genes from background
    background_genes = list(background.genes)
    if size > len(background_genes):
        raise ValueError(f"Requested size {size} exceeds background size {len(background_genes)}")
    
    sampled_genes = random.sample(background_genes, size)
    
    # Create GeneSet (using background genes as validation set)
    gene_set = GeneSet(
        gene_list=sampled_genes,
        validation_set=background.genes,
        name=f"random_{size}_{seed}",
        hgcn=True,
        format=True
    )
    
    return gene_set


def run_igea_for_library(
    gene_set: GeneSet,
    library: GeneSetLibrary,
    background: BackgroundGeneSet,
    params: Dict
) -> Optional[pd.DataFrame]:
    """
    Run iGEA for a single library and return results as DataFrame.
    
    Args:
        gene_set: Input gene set
        library: Gene set library
        background: Background gene set
        params: Analysis parameters
        
    Returns:
        DataFrame with results, or None if no results
    """
    try:
        iter_enrich = IterativeEnrichment(
            gene_set=gene_set,
            gene_set_library=library,
            background_gene_set=background,
            min_term_size=params["min_term_size"],
            max_term_size=params["max_term_size"],
            p_value_method_name=params["p_value_method"],
            p_threshold=params["p_threshold"],
            max_iterations=params["max_iterations"],
            min_overlap=params["min_overlap"],
            progress_callback=None,  # No progress callback for batch processing
        )
        
        if iter_enrich.results:
            df = iter_enrich.to_dataframe()
            # Add library name if not present
            if "Library" not in df.columns or df["Library"].isna().all():
                df["Library"] = library.name
            return df
        else:
            return None
            
    except Exception as e:
        logger.error(f"Error running iGEA for library {library.name}: {e}")
        return None


def run_single_permutation(
    args: tuple
) -> tuple[int, int, bool, Optional[str]]:
    """
    Run a single permutation (all libraries for one random gene set).
    
    Args:
        args: Tuple of (size, perm_idx, background_genes, library_paths, params, output_dir)
              Note: background_genes is a set of gene symbols (picklable)
              library_paths is a dict mapping library names to file paths
        
    Returns:
        Tuple of (size, perm_idx, success, error_message)
    """
    size, perm_idx, background_genes, library_paths, params, output_dir = args
    
    try:
        # Reload background and libraries in worker process (for multiprocessing compatibility)
        background = BackgroundGeneSet(
            str(BACKGROUNDS_DIR / "all_genes.txt"),
            name="all_genes",
            input_format="symbols"
        )
        
        libraries = {}
        for lib_name, lib_path in library_paths.items():
            library = GeneSetLibrary(str(lib_path), name=lib_name)
            # Filter terms by size
            filtered_terms = [
                t for t in library.library
                if params["min_term_size"] <= t["size"] <= params["max_term_size"]
            ]
            library.library = filtered_terms
            library.num_terms = len(filtered_terms)
            library.unique_genes = library.compute_unique_genes()
            library.size = len(library.unique_genes)
            libraries[lib_name] = library
        
        # Create random gene set
        seed = perm_idx  # Use permutation index as seed for reproducibility
        gene_set = create_random_gene_set(background, size, seed=seed)
        
        # Run iGEA for each library
        all_results = []
        for lib_name, library in libraries.items():
            df = run_igea_for_library(gene_set, library, background, params)
            if df is not None and not df.empty:
                all_results.append(df)
        
        # Combine results from all libraries
        if all_results:
            combined_df = pd.concat(all_results, ignore_index=True)
            
            # Save to TSV
            output_file = output_dir / f"permutation_{perm_idx:04d}.tsv"
            combined_df.to_csv(output_file, sep="\t", index=False)
            
            return (size, perm_idx, True, None)
        else:
            # No results from any library
            error_msg = "No results from any library"
            logger.warning(f"Size {size}, Permutation {perm_idx}: {error_msg}")
            return (size, perm_idx, False, error_msg)
            
    except Exception as e:
        error_msg = str(e)
        logger.error(f"Size {size}, Permutation {perm_idx}: Error - {error_msg}")
        return (size, perm_idx, False, error_msg)


def get_completed_permutations(output_dir: Path, size: int) -> Set[int]:
    """Get set of already completed permutation indices for a given size."""
    size_dir = output_dir / f"size_{size}"
    if not size_dir.exists():
        return set()
    
    completed = set()
    for tsv_file in size_dir.glob("permutation_*.tsv"):
        try:
            # Extract permutation index from filename
            perm_idx = int(tsv_file.stem.split("_")[1])
            completed.add(perm_idx)
        except (ValueError, IndexError):
            continue
    
    return completed


def run_permutations_for_size(
    size: int,
    n_permutations: int,
    background: BackgroundGeneSet,
    libraries: Dict[str, GeneSetLibrary],
    params: Dict,
    output_dir: Path,
    n_jobs: int = 1,
    resume: bool = True
) -> Dict[str, int]:
    """
    Run all permutations for a given gene list size.
    
    Returns:
        Dictionary with success/failure counts
    """
    size_dir = output_dir / f"size_{size}"
    size_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing gene list size: {size}")
    logger.info(f"{'='*60}")
    
    # Check for already completed permutations
    if resume:
        completed = get_completed_permutations(output_dir, size)
        remaining = set(range(1, n_permutations + 1)) - completed
        logger.info(f"Found {len(completed)} completed permutations, {len(remaining)} remaining")
    else:
        remaining = set(range(1, n_permutations + 1))
        completed = set()
    
    if not remaining:
        logger.info(f"All {n_permutations} permutations already completed for size {size}")
        return {"success": n_permutations, "failed": 0, "skipped": 0}
    
    # Prepare arguments for parallel processing
    # Convert libraries to paths for pickling (multiprocessing compatibility)
    library_paths = {
        lib_name: LIBRARIES_DIR / filename
        for lib_name, filename in LIBRARIES.items()
        if (LIBRARIES_DIR / filename).exists()
    }
    background_genes = background.genes  # Set is picklable
    
    args_list = [
        (size, perm_idx, background_genes, library_paths, params, size_dir)
        for perm_idx in sorted(remaining)
    ]
    
    # Run permutations
    stats = {"success": len(completed), "failed": 0, "skipped": 0}
    total_to_process = len(args_list)
    start_time = datetime.now()
    
    # Progress reporting frequency: every 1% or every 50 permutations, whichever is more frequent
    progress_interval = max(1, min(50, total_to_process // 100))
    
    if n_jobs == 1:
        # Sequential processing
        for i, args in enumerate(args_list, 1):
            result = run_single_permutation(args)
            _, _, success, error = result
            
            if success:
                stats["success"] += 1
            else:
                stats["failed"] += 1
            
            # Report progress at intervals or on last item
            if i % progress_interval == 0 or i == total_to_process:
                elapsed = datetime.now() - start_time
                elapsed_seconds = elapsed.total_seconds()
                rate = i / elapsed_seconds if elapsed_seconds > 0 else 0
                remaining = total_to_process - i
                eta_seconds = remaining / rate if rate > 0 else 0
                eta = datetime.now() + timedelta(seconds=eta_seconds)
                
                pct = (i / total_to_process) * 100
                logger.info(
                    f"Progress: {i}/{total_to_process} ({pct:.1f}%) | "
                    f"Elapsed: {elapsed} | "
                    f"Rate: {rate:.2f} perm/s | "
                    f"ETA: {eta.strftime('%H:%M:%S')} | "
                    f"Success: {stats['success']}, Failed: {stats['failed']}"
                )
    else:
        # Parallel processing with progress tracking
        with Pool(processes=n_jobs) as pool:
            completed_count = 0
            for i, result in enumerate(pool.imap_unordered(run_single_permutation, args_list), 1):
                size_result, perm_idx, success, error = result
                
                if success:
                    stats["success"] += 1
                else:
                    stats["failed"] += 1
                
                completed_count += 1
                
                # Report progress at intervals or on last item
                if completed_count % progress_interval == 0 or completed_count == total_to_process:
                    elapsed = datetime.now() - start_time
                    elapsed_seconds = elapsed.total_seconds()
                    rate = completed_count / elapsed_seconds if elapsed_seconds > 0 else 0
                    remaining = total_to_process - completed_count
                    eta_seconds = remaining / rate if rate > 0 else 0
                    eta = datetime.now() + timedelta(seconds=eta_seconds)
                    
                    pct = (completed_count / total_to_process) * 100
                    logger.info(
                        f"Progress: {completed_count}/{total_to_process} ({pct:.1f}%) | "
                        f"Elapsed: {elapsed} | "
                        f"Rate: {rate:.2f} perm/s | "
                        f"ETA: {eta.strftime('%H:%M:%S')} | "
                        f"Success: {stats['success']}, Failed: {stats['failed']}"
                    )
    
    size_duration = datetime.now() - start_time
    logger.info(f"Size {size} completed: {stats['success']} successful, {stats['failed']} failed in {size_duration}")
    return stats


def main():
    """Main function to run permutation distribution generation."""
    parser = argparse.ArgumentParser(
        description="Generate permutation distribution for iGEA analysis"
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=1000,
        help="Number of permutations per gene list size (default: 1000)"
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=None,
        help="Number of parallel jobs (default: number of CPU cores)"
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        default=True,
        help="Resume from already completed permutations (default: True)"
    )
    parser.add_argument(
        "--no-resume",
        action="store_false",
        dest="resume",
        help="Do not resume, overwrite existing results"
    )
    parser.add_argument(
        "--sizes",
        type=int,
        nargs="+",
        default=None,
        help="Specific gene list sizes to process (default: all sizes 50-1000)"
    )
    
    args = parser.parse_args()
    
    n_permutations = args.n_permutations
    n_jobs = args.n_jobs if args.n_jobs is not None else cpu_count()
    resume = args.resume
    sizes_to_process = args.sizes if args.sizes else GENE_LIST_SIZES
    
    logger.info("="*60)
    logger.info("iGEA Permutation Distribution Generator")
    logger.info("="*60)
    logger.info(f"Gene list sizes: {sizes_to_process}")
    logger.info(f"Permutations per size: {n_permutations}")
    logger.info(f"Parallel jobs: {n_jobs}")
    logger.info(f"Resume mode: {resume}")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"Parameters: {DEFAULT_PARAMS}")
    logger.info("="*60)
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Save configuration
    config = {
        "libraries": LIBRARIES,
        "parameters": DEFAULT_PARAMS,
        "gene_list_sizes": sizes_to_process,
        "n_permutations": n_permutations,
        "timestamp": datetime.now().isoformat()
    }
    config_file = OUTPUT_DIR / "config.json"
    with open(config_file, "w") as f:
        json.dump(config, f, indent=2)
    logger.info(f"Saved configuration to {config_file}")
    
    # Load background and libraries (once, shared across all permutations)
    logger.info("\nLoading background and libraries...")
    background = load_background()
    libraries = load_libraries()
    logger.info(f"Loaded {len(libraries)} libraries")
    
    # Process each gene list size
    all_stats = {}
    start_time = datetime.now()
    total_sizes = len(sizes_to_process)
    
    # Calculate total permutations across all sizes for overall progress
    total_permutations = 0
    for size in sizes_to_process:
        size_dir = OUTPUT_DIR / f"size_{size}"
        if resume:
            completed = get_completed_permutations(OUTPUT_DIR, size)
            remaining = n_permutations - len(completed)
        else:
            remaining = n_permutations
        total_permutations += remaining
    
    logger.info(f"\nTotal permutations to process: {total_permutations:,} across {total_sizes} sizes")
    logger.info("="*60)
    
    processed_permutations = 0
    
    for size_idx, size in enumerate(sizes_to_process, 1):
        size_start = datetime.now()
        
        # Calculate how many permutations this size will process
        size_dir = OUTPUT_DIR / f"size_{size}"
        if resume:
            completed = get_completed_permutations(OUTPUT_DIR, size)
            size_remaining = n_permutations - len(completed)
        else:
            size_remaining = n_permutations
        
        logger.info(f"\n[{size_idx}/{total_sizes}] Processing size {size} ({size_remaining} permutations)")
        logger.info("-" * 60)
        
        stats = run_permutations_for_size(
            size=size,
            n_permutations=n_permutations,
            background=background,
            libraries=libraries,
            params=DEFAULT_PARAMS,
            output_dir=OUTPUT_DIR,
            n_jobs=n_jobs,
            resume=resume
        )
        all_stats[size] = stats
        
        processed_permutations += size_remaining
        size_duration = datetime.now() - size_start
        
        # Overall progress
        overall_elapsed = datetime.now() - start_time
        overall_rate = processed_permutations / overall_elapsed.total_seconds() if overall_elapsed.total_seconds() > 0 else 0
        remaining_permutations = total_permutations - processed_permutations
        overall_eta_seconds = remaining_permutations / overall_rate if overall_rate > 0 else 0
        overall_eta = datetime.now() + timedelta(seconds=overall_eta_seconds)
        overall_pct = (processed_permutations / total_permutations) * 100 if total_permutations > 0 else 0
        
        logger.info(f"Size {size} completed in {size_duration}")
        logger.info(
            f"Overall progress: {processed_permutations:,}/{total_permutations:,} ({overall_pct:.1f}%) | "
            f"Elapsed: {overall_elapsed} | "
            f"Rate: {overall_rate:.2f} perm/s | "
            f"ETA: {overall_eta.strftime('%Y-%m-%d %H:%M:%S')}"
        )
    
    # Summary
    total_duration = datetime.now() - start_time
    logger.info("\n" + "="*60)
    logger.info("Summary")
    logger.info("="*60)
    
    total_success = sum(s["success"] for s in all_stats.values())
    total_failed = sum(s["failed"] for s in all_stats.values())
    total_skipped = sum(s.get("skipped", 0) for s in all_stats.values())
    
    logger.info(f"Total successful: {total_success}")
    logger.info(f"Total failed: {total_failed}")
    logger.info(f"Total skipped: {total_skipped}")
    logger.info(f"Total duration: {total_duration}")
    logger.info("="*60)
    
    # Save summary
    summary = {
        "total_success": total_success,
        "total_failed": total_failed,
        "total_skipped": total_skipped,
        "duration_seconds": total_duration.total_seconds(),
        "per_size_stats": all_stats
    }
    summary_file = OUTPUT_DIR / "summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Saved summary to {summary_file}")


if __name__ == "__main__":
    main()

