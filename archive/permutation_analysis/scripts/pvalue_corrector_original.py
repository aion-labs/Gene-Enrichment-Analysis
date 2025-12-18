#!/usr/bin/env python3
"""
P-value correction for iGEA based on stratified null distributions from permutations.

This module implements p-value correction that accounts for dependencies on:
- Library (GO BP, GO CC, GO MF, KEGG, Reactome)
- Iteration number
- Term size (number of genes in term)
- Overlap size (number of overlapping genes)
- Gene list size (size of input gene list)
"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple, Optional
from scipy.stats import rv_histogram
from collections import defaultdict

logger = logging.getLogger(__name__)


class iGEAPValueCorrector:
    """
    Corrects iGEA p-values using stratified empirical null distributions.
    
    The correction accounts for dependencies on library, iteration number,
    term size, overlap size, and gene list size by creating separate null
    distributions for each combination of these factors.
    
    Note: Permutation data has gene list sizes at 50, 100, 150, ..., 1000.
    For intermediate sizes, we use the nearest neighbor approach (find closest
    permutation size) or bin into ranges.
    """
    
    def __init__(
        self,
        permutation_df: pd.DataFrame,
        min_samples_per_stratum: int = 100,
        p_value_col: str = 'iteration p-value',
        overlap_col: str = 'iteration overlapping genes'
    ):
        """
        Initialize the p-value corrector.
        
        Args:
            permutation_df: DataFrame with permutation results
            min_samples_per_stratum: Minimum samples required per stratum
            p_value_col: Column name for p-values
            overlap_col: Column name for overlap information (format: "X/Y")
        """
        self.permutation_df = permutation_df.copy()
        self.min_samples = min_samples_per_stratum
        self.p_value_col = p_value_col
        self.overlap_col = overlap_col
        self.stratum_cdfs: Dict[Tuple, Optional[rv_histogram]] = {}
        self.global_cdf: Optional[rv_histogram] = None
        
        # Get available gene list sizes from permutation data
        self.available_gene_list_sizes = None
        
        # Parse permutation data
        self._parse_data()
        # Build null distributions
        self._build_null_distributions()
    
    def _parse_data(self):
        """Parse permutation data to extract features."""
        # Parse overlap information
        if self.overlap_col in self.permutation_df.columns:
            self.permutation_df['overlap_size'] = (
                self.permutation_df[self.overlap_col]
                .str.split('/')
                .str[0]
                .astype(int)
            )
            self.permutation_df['term_size'] = (
                self.permutation_df[self.overlap_col]
                .str.split('/')
                .str[1]
                .astype(int)
            )
        else:
            logger.warning(f"Column '{self.overlap_col}' not found, using defaults")
            self.permutation_df['overlap_size'] = 3
            self.permutation_df['term_size'] = 100
        
        # Parse p-values
        if self.p_value_col in self.permutation_df.columns:
            self.permutation_df['p_value'] = (
                self.permutation_df[self.p_value_col].astype(float)
            )
        else:
            raise ValueError(f"P-value column '{self.p_value_col}' not found")
        
        # Get available gene list sizes from permutation data
        if 'gene_list_size' in self.permutation_df.columns:
            self.available_gene_list_sizes = sorted(
                self.permutation_df['gene_list_size'].unique()
            )
            logger.info(f"Found gene list sizes in permutation data: {self.available_gene_list_sizes}")
        else:
            logger.warning("No 'gene_list_size' column found in permutation data")
            self.available_gene_list_sizes = []
    
    def _bin_iteration(self, iteration: int) -> int:
        """Bin iteration number."""
        if iteration == 1:
            return 1
        elif iteration <= 3:
            return 2
        elif iteration <= 6:
            return 3
        elif iteration <= 10:
            return 4
        else:
            return 5
    
    def _bin_term_size(self, term_size: int) -> str:
        """Bin term size."""
        if term_size < 50:
            return '<50'
        elif term_size < 100:
            return '50-100'
        elif term_size < 200:
            return '100-200'
        elif term_size < 300:
            return '200-300'
        else:
            return '300+'
    
    def _bin_overlap_size(self, overlap_size: int) -> str:
        """Bin overlap size."""
        if overlap_size <= 3:
            return '3'
        elif overlap_size == 4:
            return '4'
        elif overlap_size == 5:
            return '5'
        elif overlap_size <= 7:
            return '6-7'
        else:
            return '8+'
    
    def _bin_gene_list_size(self, gene_list_size: int) -> int:
        """
        Bin gene list size or find nearest available size.
        
        Permutation data has sizes: 50, 100, 150, ..., 1000.
        For intermediate sizes, we find the nearest available size.
        
        Args:
            gene_list_size: Size of the gene list
            
        Returns:
            Binned gene list size (nearest available size from permutations)
        """
        if self.available_gene_list_sizes is None or len(self.available_gene_list_sizes) == 0:
            # Fallback: bin into ranges
            if gene_list_size < 75:
                return 50
            elif gene_list_size < 125:
                return 100
            elif gene_list_size < 175:
                return 150
            elif gene_list_size < 225:
                return 200
            elif gene_list_size < 275:
                return 250
            elif gene_list_size < 325:
                return 300
            elif gene_list_size < 375:
                return 350
            elif gene_list_size < 425:
                return 400
            elif gene_list_size < 475:
                return 450
            elif gene_list_size < 525:
                return 500
            elif gene_list_size < 575:
                return 550
            elif gene_list_size < 625:
                return 600
            elif gene_list_size < 675:
                return 650
            elif gene_list_size < 725:
                return 700
            elif gene_list_size < 775:
                return 750
            elif gene_list_size < 825:
                return 800
            elif gene_list_size < 875:
                return 850
            elif gene_list_size < 925:
                return 900
            elif gene_list_size < 975:
                return 950
            else:
                return 1000
        
        # Find nearest available size
        available = np.array(self.available_gene_list_sizes)
        idx = np.abs(available - gene_list_size).argmin()
        return int(available[idx])
    
    def _create_stratum(
        self,
        library: str,
        iteration: int,
        term_size: int,
        overlap_size: int,
        gene_list_size: Optional[int] = None
    ) -> Tuple:
        """
        Create stratum identifier.
        
        Args:
            library: Library name
            iteration: Iteration number
            term_size: Number of genes in term
            overlap_size: Number of overlapping genes
            gene_list_size: Size of input gene list (optional, will be binned)
            
        Returns:
            Stratum tuple: (library, iteration_bin, term_size_bin, overlap_bin, gene_list_size_bin)
        """
        gene_list_bin = None
        if gene_list_size is not None:
            gene_list_bin = self._bin_gene_list_size(gene_list_size)
        
        return (
            library,
            self._bin_iteration(iteration),
            self._bin_term_size(term_size),
            self._bin_overlap_size(overlap_size),
            gene_list_bin
        )
    
    def _build_null_distributions(self):
        """Build empirical CDFs for each stratum."""
        df = self.permutation_df
        
        # Create strata (include gene_list_size if available)
        df['stratum'] = df.apply(
            lambda row: self._create_stratum(
                row.get('Library', 'Unknown'),
                row.get('Iteration', 1),
                row.get('term_size', 100),
                row.get('overlap_size', 3),
                row.get('gene_list_size', None)
            ),
            axis=1
        )
        
        # Build CDF for each stratum
        stratum_counts = defaultdict(int)
        for stratum, stratum_data in df.groupby('stratum'):
            p_values = stratum_data['p_value'].values
            
            if len(p_values) >= self.min_samples:
                # Create histogram (focus on range [0, 0.05] since p-values are thresholded)
                counts, bins = np.histogram(
                    p_values,
                    bins=100,
                    range=(0, 0.05),
                    density=False
                )
                
                # Avoid division by zero
                if counts.sum() > 0:
                    # Normalize to get probability density
                    counts = counts.astype(float) / counts.sum()
                    # Create CDF
                    try:
                        cdf = rv_histogram((counts, bins))
                        self.stratum_cdfs[stratum] = cdf
                        stratum_counts[stratum] = len(p_values)
                    except Exception as e:
                        logger.warning(f"Failed to create CDF for stratum {stratum}: {e}")
                        self.stratum_cdfs[stratum] = None
                else:
                    self.stratum_cdfs[stratum] = None
            else:
                # Mark as sparse
                self.stratum_cdfs[stratum] = None
                logger.debug(f"Stratum {stratum} has only {len(p_values)} samples (< {self.min_samples})")
        
        # Build global CDF as fallback
        all_p_values = df['p_value'].values
        counts, bins = np.histogram(all_p_values, bins=100, range=(0, 0.05), density=False)
        if counts.sum() > 0:
            counts = counts.astype(float) / counts.sum()
            self.global_cdf = rv_histogram((counts, bins))
        
        logger.info(f"Built {len([c for c in self.stratum_cdfs.values() if c is not None])} stratum CDFs")
        logger.info(f"Total strata: {len(self.stratum_cdfs)}")
    
    def correct_pvalue(
        self,
        p_raw: float,
        library: str,
        iteration: int,
        term_size: int,
        overlap_size: int,
        gene_list_size: Optional[int] = None
    ) -> float:
        """
        Correct a p-value using the stratified null distribution.
        
        Args:
            p_raw: Raw p-value to correct
            library: Library name
            iteration: Iteration number
            term_size: Number of genes in term
            overlap_size: Number of overlapping genes
            gene_list_size: Size of input gene list (optional but recommended)
            
        Returns:
            Corrected p-value (transformed to uniform scale under null)
        """
        stratum = self._create_stratum(library, iteration, term_size, overlap_size, gene_list_size)
        
        if stratum in self.stratum_cdfs and self.stratum_cdfs[stratum] is not None:
            # Use stratum-specific CDF
            cdf = self.stratum_cdfs[stratum]
            # Transform to uniform scale
            # If p_raw is beyond the CDF range, clamp it
            p_clamped = min(max(p_raw, 0), 0.05)
            p_corrected = cdf.cdf(p_clamped)
        else:
            # Fallback: use parent stratum or global correction
            p_corrected = self._fallback_correction(p_raw, stratum)
        
        return p_corrected
    
    def _fallback_correction(self, p_raw: float, stratum: Tuple) -> float:
        """
        Fallback correction when stratum is sparse.
        Tries parent strata (removing dimensions one at a time).
        """
        # Unpack stratum (now includes gene_list_size)
        if len(stratum) == 5:
            library, iter_bin, term_bin, overlap_bin, gene_list_bin = stratum
        else:
            # Backward compatibility
            library, iter_bin, term_bin, overlap_bin = stratum
            gene_list_bin = None
        
        # Try without gene list size bin
        if gene_list_bin is not None:
            parent_stratum = (library, iter_bin, term_bin, overlap_bin, None)
            if parent_stratum in self.stratum_cdfs and self.stratum_cdfs[parent_stratum] is not None:
                cdf = self.stratum_cdfs[parent_stratum]
                p_clamped = min(max(p_raw, 0), 0.05)
                return cdf.cdf(p_clamped)
        
        # Try without overlap bin
        parent_stratum = (library, iter_bin, term_bin, None, gene_list_bin if gene_list_bin is not None else None)
        if parent_stratum in self.stratum_cdfs and self.stratum_cdfs[parent_stratum] is not None:
            cdf = self.stratum_cdfs[parent_stratum]
            p_clamped = min(max(p_raw, 0), 0.05)
            return cdf.cdf(p_clamped)
        
        # Try without term bin
        parent_stratum = (library, iter_bin, None, None, gene_list_bin if gene_list_bin is not None else None)
        if parent_stratum in self.stratum_cdfs and self.stratum_cdfs[parent_stratum] is not None:
            cdf = self.stratum_cdfs[parent_stratum]
            p_clamped = min(max(p_raw, 0), 0.05)
            return cdf.cdf(p_clamped)
        
        # Try library-specific (any iteration, any term size, any overlap, any gene list size)
        if self.global_cdf is not None:
            p_clamped = min(max(p_raw, 0), 0.05)
            return self.global_cdf.cdf(p_clamped)
        
        # Final fallback: no correction
        logger.warning(f"No correction available for stratum {stratum}, returning raw p-value")
        return p_raw
    
    def correct_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Correct p-values in a DataFrame.
        
        Args:
            df: DataFrame with columns: Library, Iteration, 
                iteration overlapping genes, iteration p-value
                
        Returns:
            DataFrame with added 'corrected p-value' column
        """
        df = df.copy()
        
        # Parse overlap information if needed
        if 'overlap_size' not in df.columns:
            df['overlap_size'] = (
                df['iteration overlapping genes']
                .str.split('/')
                .str[0]
                .astype(int)
            )
            df['term_size'] = (
                df['iteration overlapping genes']
                .str.split('/')
                .str[1]
                .astype(int)
            )
        
        # Correct each p-value
        df['corrected p-value'] = df.apply(
            lambda row: self.correct_pvalue(
                float(row.get('iteration p-value', 1.0)),
                row.get('Library', 'Unknown'),
                int(row.get('Iteration', 1)),
                int(row.get('term_size', 100)),
                int(row.get('overlap_size', 3)),
                int(row.get('gene_list_size', None)) if 'gene_list_size' in row and pd.notna(row.get('gene_list_size')) else None
            ),
            axis=1
        )
        
        return df


def load_corrector_from_permutations(
    permutation_file: str,
    min_samples_per_stratum: int = 100
) -> iGEAPValueCorrector:
    """
    Load a corrector from permutation results file.
    
    Args:
        permutation_file: Path to merged permutation results TSV
        min_samples_per_stratum: Minimum samples per stratum
        
    Returns:
        Initialized iGEAPValueCorrector
    """
    logger.info(f"Loading permutation data from {permutation_file}")
    df = pd.read_csv(permutation_file, sep='\t')
    logger.info(f"Loaded {len(df)} permutation results")
    
    corrector = iGEAPValueCorrector(
        df,
        min_samples_per_stratum=min_samples_per_stratum
    )
    
    return corrector


if __name__ == "__main__":
    # Example usage
    import sys
    from pathlib import Path
    
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    
    # Find permutation results
    project_root = Path(__file__).parent.parent
    permutation_file = project_root / "permutation_results" / "merged_permutation_results.tsv"
    
    if not permutation_file.exists():
        print(f"Error: Permutation file not found at {permutation_file}")
        sys.exit(1)
    
    # Load corrector
    corrector = load_corrector_from_permutations(str(permutation_file))
    
    # Example: correct a p-value
    p_raw = 0.01
    library = "GO BP"
    iteration = 2
    term_size = 150
    overlap_size = 4
    gene_list_size = 73  # Example: user has 73 genes (not in permutation data)
    
    p_corrected = corrector.correct_pvalue(
        p_raw, library, iteration, term_size, overlap_size, gene_list_size
    )
    
    print(f"Raw p-value: {p_raw}")
    print(f"Gene list size: {gene_list_size} (will use nearest: {corrector._bin_gene_list_size(gene_list_size)})")
    print(f"Corrected p-value: {p_corrected:.6f}")
    print(f"Correction factor: {p_corrected / p_raw:.3f}")
