#!/usr/bin/env python3
"""
Create STRING-DB protein interaction library for enrichment analysis.

This script:
1. Downloads the STRING-DB human protein interaction file
2. Filters interactions by score >= 0.7
3. Maps Ensembl IDs to gene symbols using gene_info
4. Creates a GMT format library
5. Adds it to the enrichment system
"""

import gzip
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
STRING_DB_URL = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
PROJECT_ROOT = Path(__file__).resolve().parent
LIBRARIES_DIR = PROJECT_ROOT / "data" / "libraries"
GENE_INFO_PATH = PROJECT_ROOT / "Homo_sapiens.gene_info"
TEMP_DIR = PROJECT_ROOT / "data" / "recent_release"


def download_stringdb_file() -> Path:
    """Download the STRING-DB human protein interaction file."""
    logger.info(f"Downloading STRING-DB file from {STRING_DB_URL}")
    
    # Create temp directory
    TEMP_DIR.mkdir(exist_ok=True)
    
    # Download file
    output_path = TEMP_DIR / "9606.protein.links.v12.0.txt.gz"
    
    try:
        import requests
        response = requests.get(STRING_DB_URL, stream=True)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        logger.info(f"Downloaded {output_path}")
        return output_path
        
    except ImportError:
        logger.error("requests module not found. Please install it with: pip install requests")
        logger.info("Alternatively, you can manually download the file and place it in the temp_stringdb directory")
        raise
    except Exception as e:
        logger.error(f"Error downloading STRING-DB file: {e}")
        raise


def load_gene_info_mappings() -> Dict[str, str]:
    """
    Load gene info mappings from Homo_sapiens.gene_info.
    Maps Ensembl IDs to gene symbols.
    """
    logger.info("Loading gene info mappings...")
    
    ensembl_to_symbol = {}
    
    if not GENE_INFO_PATH.exists():
        logger.error(f"Gene info file not found at {GENE_INFO_PATH}")
        raise FileNotFoundError(f"Gene info file not found at {GENE_INFO_PATH}")
    
    try:
        with open(GENE_INFO_PATH, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                try:
                    fields = line.strip().split('\t')
                    if len(fields) < 7:
                        continue
                    
                    tax_id = fields[0]
                    if tax_id != '9606':  # Only human genes
                        continue
                    
                    entrez_id = fields[1]
                    symbol = fields[2]
                    dbxrefs = fields[5] if len(fields) > 5 else None
                    
                    # Skip entries with missing data
                    if not symbol or symbol == '-' or not dbxrefs or dbxrefs == '-':
                        continue
                    
                    # Extract Ensembl ID from dbXrefs field
                    ensembl_id = None
                    for ref in dbxrefs.split('|'):
                        if ref.startswith('Ensembl:'):
                            ensembl_id = ref.replace('Ensembl:', '')
                            break
                    
                    if ensembl_id:
                        # Store Ensembl to symbol mapping
                        ensembl_to_symbol[ensembl_id] = symbol.upper()
                    
                except Exception as e:
                    logger.warning(f"Error parsing gene_info line {line_num}: {e}")
                    continue
        
        logger.info(f"Loaded {len(ensembl_to_symbol)} Ensembl to symbol mappings")
        return ensembl_to_symbol
        
    except Exception as e:
        logger.error(f"Error loading gene_info: {e}")
        raise


def load_protein_to_symbol_mappings() -> Dict[str, str]:
    """
    Load protein to symbol mappings from STRING-DB protein info file.
    Maps Ensembl protein IDs to gene symbols.
    """
    logger.info("Loading protein to symbol mappings...")
    
    protein_to_symbol = {}
    protein_info_path = TEMP_DIR / "9606.protein.info.v12.0.txt.gz"
    
    if not protein_info_path.exists():
        logger.error(f"Protein info file not found at {protein_info_path}")
        raise FileNotFoundError(f"Protein info file not found at {protein_info_path}")
    
    try:
        with gzip.open(protein_info_path, 'rt', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                try:
                    fields = line.strip().split('\t')
                    if len(fields) < 2:
                        continue
                    
                    protein_id_full = fields[0]
                    preferred_name = fields[1]
                    
                    # Skip entries with missing data
                    if not preferred_name or preferred_name == '-':
                        continue
                    
                    # Extract protein ID (remove taxonomy prefix)
                    protein_id = protein_id_full.split('.', 1)[1] if '.' in protein_id_full else protein_id_full
                    
                    # Store protein to symbol mapping
                    protein_to_symbol[protein_id] = preferred_name.upper()
                    
                except Exception as e:
                    logger.warning(f"Error parsing protein_info line {line_num}: {e}")
                    continue
        
        logger.info(f"Loaded {len(protein_to_symbol)} protein to symbol mappings")
        return protein_to_symbol
        
    except Exception as e:
        logger.error(f"Error loading protein_info: {e}")
        raise


def parse_stringdb_interactions(file_path: Path, protein_to_symbol: Dict[str, str], min_score: float = 0.9) -> Dict[str, Set[str]]:
    """
    Parse STRING-DB interactions and create gene interaction sets.
    
    Args:
        file_path: Path to the STRING-DB file
        protein_to_symbol: Mapping from Ensembl protein IDs to gene symbols
        min_score: Minimum interaction score (default: 0.9)
    
    Returns:
        Dictionary mapping gene symbols to sets of interacting gene symbols
    """
    logger.info(f"Parsing STRING-DB interactions with min_score >= {min_score}")
    
    # Dictionary to store interactions for each gene
    gene_interactions = defaultdict(set)
    
    # Counters for statistics
    total_lines = 0
    filtered_lines = 0
    mapped_interactions = 0
    skipped_no_mapping = 0
    
    try:
        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                total_lines += 1
                
                if line.startswith('protein1'):
                    continue  # Skip header
                
                try:
                    fields = line.strip().split()
                    if len(fields) < 3:
                        continue
                    
                    protein1_full = fields[0]
                    protein2_full = fields[1]
                    score = float(fields[2]) / 1000.0  # Convert from integer to float
                    
                    # Filter by score
                    if score < min_score:
                        continue
                    
                    filtered_lines += 1
                    
                    # Extract Ensembl protein IDs (remove taxonomy prefix)
                    protein1 = protein1_full.split('.', 1)[1] if '.' in protein1_full else protein1_full
                    protein2 = protein2_full.split('.', 1)[1] if '.' in protein2_full else protein2_full
                    
                    # Map protein IDs to symbols
                    symbol1 = protein_to_symbol.get(protein1)
                    symbol2 = protein_to_symbol.get(protein2)
                    
                    if symbol1 and symbol2:
                        # Add interaction (bidirectional)
                        gene_interactions[symbol1].add(symbol2)
                        gene_interactions[symbol2].add(symbol1)
                        mapped_interactions += 1
                    else:
                        skipped_no_mapping += 1
                    
                    # Progress indicator
                    if line_num % 1000000 == 0:
                        logger.info(f"Processed {line_num:,} lines...")
                
                except Exception as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue
        
        logger.info(f"Parsing complete:")
        logger.info(f"  Total lines: {total_lines:,}")
        logger.info(f"  Filtered lines (score >= {min_score}): {filtered_lines:,}")
        logger.info(f"  Mapped interactions: {mapped_interactions:,}")
        logger.info(f"  Skipped (no mapping): {skipped_no_mapping:,}")
        logger.info(f"  Unique genes with interactions: {len(gene_interactions):,}")
        
        return dict(gene_interactions)
        
    except Exception as e:
        logger.error(f"Error parsing STRING-DB file: {e}")
        raise


def create_gmt_library(gene_interactions: Dict[str, Set[str]], output_path: Path) -> None:
    """
    Create a GMT format library from gene interactions.
    
    Args:
        gene_interactions: Dictionary mapping gene symbols to sets of interacting genes
        output_path: Path to output GMT file
    """
    logger.info(f"Creating GMT library at {output_path}")
    
    # Create libraries directory if it doesn't exist
    LIBRARIES_DIR.mkdir(parents=True, exist_ok=True)
    
    # Sort genes for consistent output
    sorted_genes = sorted(gene_interactions.keys())
    
    with open(output_path, 'w', encoding='utf-8') as f:
        for gene in sorted_genes:
            interactors = gene_interactions[gene]
            
            # Remove self-interaction
            interactors.discard(gene)
            
            if interactors:
                # Sort interactors for consistent output
                sorted_interactors = sorted(interactors)
                
                # Write GMT line: name, description, genes...
                line = f"{gene}\tProtein Interaction for {gene}\t" + "\t".join(sorted_interactors)
                f.write(line + "\n")
    
    logger.info(f"Created GMT library with {len(sorted_genes)} gene sets")


def update_alias_json(gmt_filename: str) -> None:
    """Update the alias.json file to include the new STRING-DB library."""
    alias_file = LIBRARIES_DIR / "alias.json"
    
    if not alias_file.exists():
        logger.warning("alias.json not found, creating new one")
        aliases = []
    else:
        # Read existing aliases
        with open(alias_file, 'r') as f:
            aliases = json.load(f)
    
    # Check if the library already exists
    existing_files = {alias['file'] for alias in aliases}
    
    if gmt_filename not in existing_files:
        new_alias = {
            "name": "Protein Interaction",
            "file": gmt_filename,
            "active": False  # Default to inactive
        }
        aliases.append(new_alias)
        logger.info(f"Added new alias: {new_alias['name']} -> {gmt_filename}")
    else:
        logger.info(f"Library {gmt_filename} already exists in alias.json")
    
    # Write updated aliases
    with open(alias_file, 'w') as f:
        json.dump(aliases, f, indent=4)
    
    logger.info(f"Updated alias.json with {len(aliases)} libraries")


def cleanup_temp_files():
    """Clean up temporary files."""
    if TEMP_DIR.exists():
        import shutil
        shutil.rmtree(TEMP_DIR)
        logger.info(f"Cleaned up temporary directory {TEMP_DIR}")


def main():
    """Main function to create STRING-DB library."""
    try:
        logger.info("Starting STRING-DB library creation...")
        
        # Step 1: Download STRING-DB files (if not already present)
        stringdb_file = TEMP_DIR / "9606.protein.links.v12.0.txt.gz"
        protein_info_file = TEMP_DIR / "9606.protein.info.v12.0.txt.gz"
        
        if not stringdb_file.exists():
            logger.info("STRING-DB protein links file not found, downloading...")
            download_stringdb_file()
        
        if not protein_info_file.exists():
            logger.info("STRING-DB protein info file not found, downloading...")
            # Download protein info file
            try:
                import requests
                response = requests.get("https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz", stream=True)
                response.raise_for_status()
                
                with open(protein_info_file, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                
                logger.info(f"Downloaded {protein_info_file}")
            except ImportError:
                logger.error("requests module not found. Please install it with: pip install requests")
                logger.info("Alternatively, you can manually download the protein info file and place it in the temp_stringdb directory")
                raise
            except Exception as e:
                logger.error(f"Error downloading protein info file: {e}")
                raise
        
        # Step 2: Load protein to symbol mappings
        protein_to_symbol = load_protein_to_symbol_mappings()
        
        # Step 3: Parse interactions
        gene_interactions = parse_stringdb_interactions(stringdb_file, protein_to_symbol, min_score=0.9)
        
        # Step 4: Create GMT library
        gmt_filename = "stringdb_protein_interactions.gmt"
        gmt_path = LIBRARIES_DIR / gmt_filename
        create_gmt_library(gene_interactions, gmt_path)
        
        # Step 5: Update alias.json
        update_alias_json(gmt_filename)
        
        # Step 6: Cleanup (optional - comment out if you want to keep temp files)
        # cleanup_temp_files()
        
        logger.info("STRING-DB library creation completed successfully!")
        
    except Exception as e:
        logger.error(f"Error creating STRING-DB library: {e}")
        raise


if __name__ == "__main__":
    import json
    main()
