from pathlib import Path
from typing import Set, Optional

from gene_converter import GeneConverter


class BackgroundGeneSet:
    """
    A class to store a list of genes and their size.
    """

    def __init__(
        self, 
        background_file_path: str, 
        name: str = "", 
        organism: str = "Homo Sapiens",
        input_format: str = "symbols",
        skip_validation: bool = False
    ) -> None:
        """
        Initialize BackgroundGeneList object with a list of genes.

        Args:
            background_file_path: Path to the background file
            name: Name for the background gene list
            organism: Organism name
            input_format: Either 'symbols' or 'entrez_ids'
            skip_validation: If True, skip gene symbol validation (faster, use for trusted sources)
        """
        self.genes: Set[str] = self._load_from_file(background_file_path, input_format, skip_validation)
        self.size: int = len(self.genes)
        self.name = name if name else Path(background_file_path).stem
        self.organism = organism
        self.input_format = input_format

    def _load_from_file(self, background_file_path: str, input_format: str = "symbols", skip_validation: bool = False) -> Set[str]:
        """
        Load background genes from a file with optional Entrez ID conversion
        
        Args:
            background_file_path: Path to the background file
            input_format: Either 'symbols' or 'entrez_ids'
            skip_validation: If True, skip gene symbol validation (faster, use for trusted sources)

        Returns:
            Set of gene symbols representing the background
        """
        # Read raw lines from file
        with open(background_file_path, "r") as f:
            raw_lines = [line.strip() for line in f.readlines() if line.strip()]
        
        if input_format == "entrez_ids":
            # Convert Entrez IDs to symbols
            converter = GeneConverter()
            converted_symbols = []
            unrecognized_entrez = []
            
            for line in raw_lines:
                gene_id = line.strip()
                if not gene_id:
                    continue
                
                symbol = converter.get_symbol(gene_id)
                if symbol:
                    converted_symbols.append(symbol)
                else:
                    unrecognized_entrez.append(gene_id)
            
            # Log conversion results
            if unrecognized_entrez:
                print(f"Warning: {len(unrecognized_entrez)} Entrez IDs not found in database: {unrecognized_entrez[:10]}{'...' if len(unrecognized_entrez) > 10 else ''}")
            
            return set(converted_symbols)
        else:
            # Skip validation if requested (faster, use for trusted sources like permutation tests)
            if skip_validation:
                import logging
                logger = logging.getLogger(__name__)
                logger.debug(f"Skipping gene validation for background file (skip_validation=True)")
                # Return all non-empty lines as valid symbols (uppercase for consistency)
                return set(line.strip().upper() for line in raw_lines if line.strip())
            
            # Validate symbols
            converter = GeneConverter()
            
            # Check if converter has loaded gene data
            # If gene info file wasn't loaded, converter will be empty and all symbols will be rejected
            stats = converter.get_stats()
            has_gene_data = stats.get('symbol_mappings', 0) > 0
            
            if not has_gene_data:
                # Gene converter doesn't have data - likely gene info file not loaded
                # For background files, we'll accept all symbols as-is (trusted source)
                # This is safer than rejecting everything
                import logging
                logger = logging.getLogger(__name__)
                logger.warning(
                    f"GeneConverter has no loaded data (gene info file may not be loaded). "
                    f"Accepting all {len(raw_lines)} symbols from background file without validation."
                )
                # Return all non-empty lines as valid symbols
                return set(line.strip().upper() for line in raw_lines if line.strip())
            
            valid_symbols = []
            invalid_symbols = []
            
            for line in raw_lines:
                gene_id = line.strip()
                if not gene_id:
                    continue
                
                if converter.is_symbol(gene_id):
                    valid_symbols.append(gene_id)
                else:
                    invalid_symbols.append(gene_id)
            
            # Log validation results
            if invalid_symbols:
                import logging
                logger = logging.getLogger(__name__)
                logger.warning(
                    f"{len(invalid_symbols)} symbols not found in database: {invalid_symbols[:10]}{'...' if len(invalid_symbols) > 10 else ''}"
                )
            
            # If all symbols were invalid, this might indicate a problem
            if len(valid_symbols) == 0 and len(raw_lines) > 0:
                import logging
                logger = logging.getLogger(__name__)
                logger.error(
                    f"All {len(raw_lines)} symbols from background file were rejected by GeneConverter. "
                    f"This likely indicates the gene info file is not loaded correctly. "
                    f"First 5 symbols: {[line.strip() for line in raw_lines[:5]]}"
                )
            
            return set(valid_symbols)

    def has_gene(self, gene: str) -> bool:
        """
        Check if the given gene is present in the BackgroundGenes.

        Args:
            gene: A gene name.

        Returns:
            True if the gene is present, False otherwise.
        """
        return gene in self.genes
