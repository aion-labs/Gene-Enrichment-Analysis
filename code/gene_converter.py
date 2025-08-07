import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class GeneConverter:
    """
    Converter between Entrez IDs and gene symbols using NCBI gene info data.
    """
    
    def __init__(self, gene_info_path: Optional[str] = None):
        """
        Initialize the gene converter.
        
        Args:
            gene_info_path: Path to the Homo_sapiens.gene_info file. 
                          If None, will look for it in the project root.
        """
        if gene_info_path is None:
            # Look for the file in the project root
            project_root = Path(__file__).resolve().parent.parent
            gene_info_path = project_root / "Homo_sapiens.gene_info"
        
        self.gene_info_path = Path(gene_info_path)
        self.entrez_to_symbol: Dict[str, str] = {}
        self.symbol_to_entrez: Dict[str, str] = {}
        self.synonyms_to_symbol: Dict[str, str] = {}
        
        if self.gene_info_path.exists():
            self._load_gene_info()
        else:
            logger.warning(f"Gene info file not found at {self.gene_info_path}")
    
    def _load_gene_info(self) -> None:
        """Load gene information from the NCBI gene info file."""
        logger.info(f"Loading gene info from {self.gene_info_path}")
        
        with open(self.gene_info_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                try:
                    fields = line.strip().split('\t')
                    if len(fields) < 3:
                        continue
                    
                    tax_id = fields[0]
                    if tax_id != '9606':  # Only human genes
                        continue
                    
                    entrez_id = fields[1]
                    symbol = fields[2]
                    
                    # Skip entries with missing or invalid data
                    if not entrez_id or not symbol or symbol == '-':
                        continue
                    
                    # Store primary mappings
                    self.entrez_to_symbol[entrez_id] = symbol
                    self.symbol_to_entrez[symbol] = entrez_id
                    
                    # Store synonyms if available
                    if len(fields) > 4 and fields[4] != '-':
                        synonyms = fields[4].split('|')
                        for synonym in synonyms:
                            if synonym and synonym != '-':
                                self.synonyms_to_symbol[synonym] = symbol
                
                except Exception as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue
        
        logger.info(f"Loaded {len(self.entrez_to_symbol)} human Entrez ID mappings (taxonomy ID 9606)")
        logger.info(f"Loaded {len(self.synonyms_to_symbol)} synonym mappings")
        
        # Debug: Show a few examples
        if self.entrez_to_symbol:
            sample_entrez = list(self.entrez_to_symbol.keys())[:3]
            logger.info(f"Sample Entrez IDs: {sample_entrez}")
            logger.info(f"Sample symbols: {[self.entrez_to_symbol[e] for e in sample_entrez]}")
    
    def convert_input(self, input_text: str) -> Tuple[List[str], List[str], List[str]]:
        """
        Convert input text to gene symbols, handling both Entrez IDs and gene symbols.
        
        Args:
            input_text: Text containing gene identifiers (one per line)
            
        Returns:
            Tuple of (converted_symbols, unrecognized_entrez, unrecognized_symbols)
        """
        lines = [line.strip() for line in input_text.split('\n') if line.strip()]
        
        converted_symbols = []
        unrecognized_entrez = []
        unrecognized_symbols = []
        
        for line in lines:
            # Remove any extra whitespace and common separators
            gene_id = re.sub(r'[,\s]+', '', line)
            
            if not gene_id:
                continue
            
            # Try to convert
            symbol = self.entrez_to_symbol.get(gene_id)
            if symbol:
                # It's an Entrez ID
                converted_symbols.append(symbol)
            elif gene_id in self.symbol_to_entrez:
                # It's already a valid symbol
                converted_symbols.append(gene_id)
            elif gene_id in self.synonyms_to_symbol:
                # It's a synonym
                converted_symbols.append(self.synonyms_to_symbol[gene_id])
            else:
                # Check if it looks like an Entrez ID (numeric)
                if gene_id.isdigit():
                    unrecognized_entrez.append(gene_id)
                else:
                    unrecognized_symbols.append(gene_id)
        
        return converted_symbols, unrecognized_entrez, unrecognized_symbols
    
    def get_symbol(self, entrez_id: str) -> Optional[str]:
        """Get gene symbol for an Entrez ID."""
        return self.entrez_to_symbol.get(entrez_id)
    
    def get_entrez(self, symbol: str) -> Optional[str]:
        """Get Entrez ID for a gene symbol."""
        return self.symbol_to_entrez.get(symbol)
    
    def is_entrez_id(self, gene_id: str) -> bool:
        """Check if a string looks like an Entrez ID."""
        return gene_id.isdigit() and gene_id in self.entrez_to_symbol
    
    def is_symbol(self, gene_id: str) -> bool:
        """Check if a string is a valid gene symbol."""
        return gene_id in self.symbol_to_entrez or gene_id in self.synonyms_to_symbol
    
    def get_stats(self) -> Dict[str, int]:
        """Get statistics about the loaded mappings."""
        return {
            'entrez_mappings': len(self.entrez_to_symbol),
            'symbol_mappings': len(self.symbol_to_entrez),
            'synonym_mappings': len(self.synonyms_to_symbol)
        } 