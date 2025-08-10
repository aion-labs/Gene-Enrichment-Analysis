#!/usr/bin/env python3
"""
MSigDB Update Script

This script updates the MSigDB gene set libraries in the data/libraries directory.
It can handle both zip file input and direct downloads from MSigDB.

Usage:
    python update_msigdb.py --zip-file path/to/msigdb_v2025.1.Hs_files_to_download_locally.zip
    python update_msigdb.py --download (requires login credentials)
"""

import argparse
import json
import logging
import os
import shutil
import zipfile
import gzip
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Set, Dict

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
MSIGDB_BASE_URL = "https://www.gsea-msigdb.org"
MSIGDB_DOWNLOAD_URL = "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/msigdb_v2025.1.Hs_files_to_download_locally.zip"
NCBI_GENE_INFO_URL = "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
NCBI_GENE_HISTORY_URL = "https://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz"
PROJECT_ROOT = Path(__file__).resolve().parent
LIBRARIES_DIR = PROJECT_ROOT / "data" / "libraries"
BACKUP_DIR = PROJECT_ROOT / "data" / "libraries" / "backup"
BACKGROUNDS_DIR = PROJECT_ROOT / "data" / "backgrounds"
RECENT_RELEASE_DIR = PROJECT_ROOT / "data" / "recent_release"


def download_ncbi_files() -> tuple[Optional[Path], Optional[Path]]:
    """Download Homo_sapiens.gene_info.gz and gene_history.gz from NCBI."""
    logger.info("Downloading NCBI gene information files...")
    
    try:
        import requests
    except ImportError:
        logger.error("requests module not found. Please install it with: pip install requests")
        return None, None
    
    # Create recent_release directory if it doesn't exist
    RECENT_RELEASE_DIR.mkdir(parents=True, exist_ok=True)
    
    gene_info_path = RECENT_RELEASE_DIR / "Homo_sapiens.gene_info.gz"
    gene_history_path = RECENT_RELEASE_DIR / "gene_history.gz"
    
    # Download Homo_sapiens.gene_info.gz
    logger.info(f"Downloading {NCBI_GENE_INFO_URL}...")
    try:
        response = requests.get(NCBI_GENE_INFO_URL, stream=True)
        if response.status_code == 200:
            with open(gene_info_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            logger.info(f"Downloaded {gene_info_path}")
        else:
            logger.error(f"Failed to download gene_info: {response.status_code}")
            return None, None
    except Exception as e:
        logger.error(f"Error downloading gene_info: {e}")
        return None, None
    
    # Download gene_history.gz
    logger.info(f"Downloading {NCBI_GENE_HISTORY_URL}...")
    try:
        response = requests.get(NCBI_GENE_HISTORY_URL, stream=True)
        if response.status_code == 200:
            with open(gene_history_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            logger.info(f"Downloaded {gene_history_path}")
        else:
            logger.error(f"Failed to download gene_history: {response.status_code}")
            return gene_info_path, None
    except Exception as e:
        logger.error(f"Error downloading gene_history: {e}")
        return gene_info_path, None
    
    return gene_info_path, gene_history_path


def extract_gene_info_files(gene_info_gz_path: Optional[Path], gene_history_gz_path: Optional[Path]) -> tuple[Optional[Path], Optional[Path]]:
    """Extract the gzipped NCBI files."""
    logger.info("Extracting NCBI gene information files...")
    
    gene_info_path = None
    gene_history_path = None
    
    # Extract gene_info
    if gene_info_gz_path and gene_info_gz_path.exists():
        gene_info_path = RECENT_RELEASE_DIR / "Homo_sapiens.gene_info"
        try:
            with gzip.open(gene_info_gz_path, 'rt', encoding='utf-8') as f_in:
                with open(gene_info_path, 'w', encoding='utf-8') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info(f"Extracted {gene_info_path}")
        except Exception as e:
            logger.error(f"Error extracting gene_info: {e}")
            gene_info_path = None
    
    # Extract gene_history and filter to human only
    if gene_history_gz_path and gene_history_gz_path.exists():
        gene_history_path = RECENT_RELEASE_DIR / "gene_history"
        try:
            logger.info("Extracting gene_history and filtering to human taxonomy (9606)...")
            human_count = 0
            total_count = 0
            
            with gzip.open(gene_history_gz_path, 'rt', encoding='utf-8') as f_in:
                with open(gene_history_path, 'w', encoding='utf-8') as f_out:
                    for line_num, line in enumerate(f_in, 1):
                        total_count += 1
                        
                        # Write header lines
                        if line.startswith('#'):
                            f_out.write(line)
                            continue
                        
                        # Check if it's human taxonomy (9606)
                        try:
                            fields = line.strip().split('\t')
                            if len(fields) >= 1 and fields[0] == '9606':
                                f_out.write(line)
                                human_count += 1
                        except Exception as e:
                            logger.warning(f"Error parsing gene_history line {line_num}: {e}")
                            continue
                        
                        # Progress indicator for large files
                        if line_num % 1000000 == 0:
                            logger.info(f"Processed {line_num:,} lines from gene_history...")
            
            logger.info(f"Extracted and filtered {gene_history_path}")
            logger.info(f"Retained {human_count:,} human entries out of {total_count:,} total entries")
            
        except Exception as e:
            logger.error(f"Error extracting gene_history: {e}")
            gene_history_path = None
    
    return gene_info_path, gene_history_path


def load_gene_history(gene_history_path: Path) -> Dict[str, str]:
    """Load gene history and create mapping from old symbols to current symbols."""
    logger.info("Loading gene history...")
    
    old_to_current = {}
    
    try:
        with open(gene_history_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                try:
                    fields = line.strip().split('\t')
                    if len(fields) < 4:
                        continue
                    
                    tax_id = fields[0]
                    if tax_id != '9606':  # Only human genes
                        continue
                    
                    old_symbol = fields[2]
                    current_symbol = fields[3]
                    
                    # Skip entries with missing or invalid data
                    if not old_symbol or not current_symbol or old_symbol == '-' or current_symbol == '-':
                        continue
                    
                    # Only add if there's a unique mapping (no conflicts)
                    if old_symbol not in old_to_current:
                        old_to_current[old_symbol] = current_symbol
                    elif old_to_current[old_symbol] != current_symbol:
                        # Conflict found - remove this mapping
                        del old_to_current[old_symbol]
                        logger.debug(f"Removed conflicting mapping for {old_symbol}")
                
                except Exception as e:
                    logger.warning(f"Error parsing gene_history line {line_num}: {e}")
                    continue
        
        logger.info(f"Loaded {len(old_to_current)} unique old-to-current symbol mappings")
        return old_to_current
        
    except Exception as e:
        logger.error(f"Error loading gene_history: {e}")
        return {}


def validate_and_map_genes(genes: Set[str], gene_info_path: Path, gene_history_path: Path) -> tuple[Set[str], Set[str]]:
    """
    Validate genes against Homo_sapiens.gene_info and use gene_history for mapping.
    
    Returns:
        Tuple of (valid_genes, invalid_genes)
    """
    logger.info("Validating and mapping genes...")
    
    # Load current gene symbols from gene_info
    current_symbols = set()
    try:
        with open(gene_info_path, 'r', encoding='utf-8') as f:
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
                    
                    symbol = fields[2]
                    if symbol and symbol != '-':
                        current_symbols.add(symbol.upper())
                
                except Exception as e:
                    logger.warning(f"Error parsing gene_info line {line_num}: {e}")
                    continue
        
        logger.info(f"Loaded {len(current_symbols)} current gene symbols")
        
    except Exception as e:
        logger.error(f"Error loading gene_info: {e}")
        return set(), genes
    
    # Load gene history mappings
    old_to_current = load_gene_history(gene_history_path) if gene_history_path else {}
    
    valid_genes = set()
    invalid_genes = set()
    mapped_count = 0
    
    for gene in genes:
        gene_upper = gene.upper()
        
        # First try direct match
        if gene_upper in current_symbols:
            valid_genes.add(gene_upper)
        else:
            # Try mapping through gene history
            if gene_upper in old_to_current:
                mapped_symbol = old_to_current[gene_upper].upper()
                if mapped_symbol in current_symbols:
                    valid_genes.add(mapped_symbol)
                    mapped_count += 1
                    logger.debug(f"Mapped {gene_upper} -> {mapped_symbol}")
                else:
                    invalid_genes.add(gene_upper)
            else:
                invalid_genes.add(gene_upper)
    
    logger.info(f"Validated {len(valid_genes)} genes as valid")
    logger.info(f"Mapped {mapped_count} genes using gene history")
    if invalid_genes:
        logger.warning(f"Found {len(invalid_genes)} invalid genes")
        logger.debug(f"Sample invalid genes: {list(invalid_genes)[:10]}")
    
    return valid_genes, invalid_genes


def create_backup() -> None:
    """Create a backup of the current libraries before updating."""
    if not BACKUP_DIR.exists():
        BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = BACKUP_DIR / f"backup_{timestamp}"
    
    if LIBRARIES_DIR.exists():
        shutil.copytree(LIBRARIES_DIR, backup_path, ignore=shutil.ignore_patterns('backup'))
        logger.info(f"Created backup at {backup_path}")
    else:
        logger.warning("No existing libraries directory found to backup")


def extract_zip_file(zip_path: Path) -> Path:
    """Extract the MSigDB zip file to a temporary directory."""
    temp_dir = PROJECT_ROOT / "temp_msigdb_extract"
    
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(exist_ok=True)
    
    logger.info(f"Extracting {zip_path} to {temp_dir}")
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)
    
    return temp_dir


def find_gmt_files(extract_dir: Path) -> List[Path]:
    """Find C2 and C5 .gmt files in the extracted directory (excluding c2.all and c5.all)."""
    all_gmt_files = list(extract_dir.rglob("*.gmt"))
    
    # Filter for C2 and C5 gene sets only (symbols files), excluding c2.all and c5.all
    c2_c5_files = []
    for gmt_file in all_gmt_files:
        filename = gmt_file.name.lower()
        # Check if it's a human symbols file and starts with c2 or c5, but exclude c2.all and c5.all
        if ("hs.symbols.gmt" in filename and 
            (filename.startswith("c2") or filename.startswith("c5")) and
            not filename.startswith("c2.all") and 
            not filename.startswith("c5.all")):
            c2_c5_files.append(gmt_file)
    
    logger.info(f"Found {len(all_gmt_files)} total .gmt files")
    logger.info(f"Filtered to {len(c2_c5_files)} C2 and C5 .gmt files (excluding c2.all and c5.all)")
    
    # Log the specific files being processed
    for gmt_file in c2_c5_files:
        logger.debug(f"Will process: {gmt_file.name}")
    
    return c2_c5_files


def collect_all_genes_from_gmt_files(gmt_files: List[Path]) -> Set[str]:
    """Collect all unique genes from C2 and C5 GMT files (excluding c2.all and c5.all)."""
    all_genes = set()
    
    logger.info("Collecting all unique genes from C2 and C5 GMT files (excluding c2.all and c5.all)...")
    
    for gmt_file in gmt_files:
        if "Hs.symbols.gmt" in gmt_file.name:
            logger.debug(f"Processing {gmt_file.name}")
            try:
                with open(gmt_file, 'r', encoding='utf-8') as f:
                    for line_num, line in enumerate(f, 1):
                        line = line.strip()
                        if not line:
                            continue
                        
                        # Split by tab and skip the first two columns (name and description)
                        parts = line.split('\t')
                        if len(parts) < 3:
                            continue
                        
                        # Add all genes from this line
                        genes = parts[2:]
                        for gene in genes:
                            if gene and gene.strip():
                                all_genes.add(gene.strip().upper())
                        
                        # Progress indicator for large files
                        if line_num % 10000 == 0:
                            logger.debug(f"Processed {line_num} lines from {gmt_file.name}")
                            
            except Exception as e:
                logger.warning(f"Error processing {gmt_file.name}: {e}")
    
    logger.info(f"Collected {len(all_genes)} unique genes from C2 and C5 GMT files (excluding c2.all and c5.all)")
    return all_genes


def update_background_gene_lists(valid_genes: Set[str]) -> None:
    """Update background gene lists with the validated genes from C2 and C5 gene sets (excluding c2.all and c5.all)."""
    logger.info("Updating background gene lists...")
    
    # Create backgrounds directory if it doesn't exist
    BACKGROUNDS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Sort genes for consistent output
    sorted_genes = sorted(valid_genes)
    
    # Create the main background file with all genes
    background_file = BACKGROUNDS_DIR / "all_genes.txt"
    with open(background_file, 'w') as f:
        for gene in sorted_genes:
            f.write(f"{gene}\n")
    
    logger.info(f"Created background file: {background_file} with {len(sorted_genes)} genes")
    
    # Update alias.json for backgrounds
    update_background_alias_json(background_file.name)


def update_background_alias_json(background_filename: str) -> None:
    """Update the background alias.json file."""
    alias_file = BACKGROUNDS_DIR / "alias.json"
    
    if not alias_file.exists():
        logger.warning("background alias.json not found, creating new one")
        aliases = []
    else:
        # Read existing aliases
        with open(alias_file, 'r') as f:
            aliases = json.load(f)
    
    # Clear existing aliases and create new one
    aliases = []
    
    new_alias = {
        "name": "All genes",
        "file": background_filename,
        "active": True
    }
    aliases.append(new_alias)
    logger.info(f"Added background alias: {new_alias['name']} -> {background_filename}")
    
    # Write updated aliases
    with open(alias_file, 'w') as f:
        json.dump(aliases, f, indent=4)
    
    logger.info(f"Updated background alias.json with {len(aliases)} backgrounds")


def update_libraries(gmt_files: List[Path]) -> None:
    """Update the libraries directory with new C2 and C5 .gmt files (excluding c2.all and c5.all)."""
    # Create libraries directory if it doesn't exist
    LIBRARIES_DIR.mkdir(parents=True, exist_ok=True)
    
    updated_files = []
    skipped_files = []
    
    for gmt_file in gmt_files:
        # Get the filename
        filename = gmt_file.name
        
        # Check if it's a human symbols file (Hs.symbols.gmt) and C2 or C5, but exclude c2.all and c5.all
        if ("Hs.symbols.gmt" in filename and 
            ("c2" in filename.lower() or "c5" in filename.lower()) and
            not filename.lower().startswith("c2.all") and 
            not filename.lower().startswith("c5.all")):
            target_path = LIBRARIES_DIR / filename
            
            # Copy the file
            shutil.copy2(gmt_file, target_path)
            updated_files.append(filename)
            logger.info(f"Updated {filename}")
        else:
            skipped_files.append(filename)
            logger.debug(f"Skipped {filename} (not a C2/C5 human symbols file or is c2.all/c5.all)")
    
    logger.info(f"Updated {len(updated_files)} C2 and C5 files (excluding c2.all and c5.all)")
    if skipped_files:
        logger.info(f"Skipped {len(skipped_files)} non-C2/C5 files or c2.all/c5.all files")


def generate_library_name(filename: str) -> str:
    """Generate a human-readable name from a filename."""
    # Remove .gmt extension
    name = filename.replace('.gmt', '')
    
    # Common patterns
    patterns = {
        'h.all': 'H: Hallmark gene sets',
        'c1.all': 'C1: Positional gene sets',
        'c2.all': 'C2: Curated gene sets',
        'c2.cgp': 'C2: CGP: Chemical & genetic perturbations',
        'c2.cp': 'C2: CP: Canonical pathways',
        'c2.cp.biocarta': 'C2: CP: BioCarta',
        'c2.cp.kegg_medicus': 'C2: CP: KEGG MEDICUS',
        'c2.cp.kegg_legacy': 'C2: CP: KEGG Legacy',
        'c2.cp.pid': 'C2: CP: PID',
        'c2.cp.reactome': 'C2: CP: Reactome',
        'c2.cp.wikipathways': 'C2: CP: WikiPathways',
        'c3.all': 'C3: Regulatory target gene sets',
        'c3.mir': 'C3: MIR: microRNA targets',
        'c3.mir.mirdb': 'C3: MIR: miRDB',
        'c3.mir.mir_legacy': 'C3: MIR: Legacy',
        'c3.tft': 'C3: TFT: Transcription factor targets',
        'c3.tft.gtrd': 'C3: TFT: GTRD',
        'c3.tft.tft_legacy': 'C3: TFT: Legacy',
        'c4.all': 'C4: Computational gene sets',
        'c4.3ca': 'C4: 3CA: Cancer cell lines',
        'c4.cgn': 'C4: CGN: Cancer gene neighborhoods',
        'c4.cm': 'C4: CM: Cancer modules',
        'c5.all': 'C5: GO gene sets',
        'c5.go': 'C5: GO: Gene Ontology',
        'c5.go.bp': 'C5: GO: Biological Process',
        'c5.go.cc': 'C5: GO: Cellular Component',
        'c5.go.mf': 'C5: GO: Molecular Function',
        'c5.hpo': 'C5: HPO: Human Phenotype Ontology',
        'c6.all': 'C6: Oncogenic signatures',
        'c7.all': 'C7: Immunologic signatures',
        'c7.immunesigdb': 'C7: ImmuneSigDB',
        'c7.vax': 'C7: VAX: Vaccine response',
        'c8.all': 'C8: Cell type signatures',
        'msigdb': 'All MSigDB gene sets'
    }
    
    for pattern, display_name in patterns.items():
        if pattern in name.lower():
            return display_name
    
    # If no pattern matches, create a generic name
    return f"Library: {name}"


def update_alias_json() -> None:
    """Update the alias.json file with new library information."""
    alias_file = LIBRARIES_DIR / "alias.json"
    
    if not alias_file.exists():
        logger.warning("alias.json not found, creating new one")
        aliases = []
    else:
        # Read existing aliases
        with open(alias_file, 'r') as f:
            aliases = json.load(f)
    
    # Get all .gmt files in the libraries directory
    gmt_files = list(LIBRARIES_DIR.glob("*.gmt"))
    
    # Create a mapping of existing aliases
    existing_files = {alias['file'] for alias in aliases}
    
    # Add new files
    for gmt_file in gmt_files:
        if gmt_file.name not in existing_files:
            # Generate a human-readable name
            name = generate_library_name(gmt_file.name)
            
            new_alias = {
                "name": name,
                "file": gmt_file.name,
                "active": False  # Default to inactive
            }
            aliases.append(new_alias)
            logger.info(f"Added new alias: {name} -> {gmt_file.name}")
    
    # Write updated aliases
    with open(alias_file, 'w') as f:
        json.dump(aliases, f, indent=4)
    
    logger.info(f"Updated alias.json with {len(aliases)} libraries")


def download_msigdb(email: str, password: str) -> Optional[Path]:
    """Download MSigDB files directly (requires login)."""
    try:
        import requests
    except ImportError:
        logger.error("requests module not found. Please install it with: pip install requests")
        return None
    
    session = requests.Session()
    
    # First, we need to login
    login_url = "https://www.gsea-msigdb.org/gsea/login.jsp"
    login_data = {
        'email': email,
        'password': password
    }
    
    try:
        response = session.post(login_url, data=login_data)
        if response.status_code != 200:
            logger.error("Failed to login to MSigDB")
            return None
        
        # Now download the file
        logger.info("Downloading MSigDB files...")
        download_response = session.get(MSIGDB_DOWNLOAD_URL, stream=True)
        
        if download_response.status_code == 200:
            zip_path = PROJECT_ROOT / "msigdb_download.zip"
            with open(zip_path, 'wb') as f:
                for chunk in download_response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            logger.info(f"Downloaded {zip_path}")
            return zip_path
        else:
            logger.error(f"Failed to download: {download_response.status_code}")
            return None
            
    except Exception as e:
        logger.error(f"Error during download: {e}")
        return None


def cleanup_temp_files(temp_dir: Optional[Path] = None, zip_path: Optional[Path] = None):
    """Clean up temporary files."""
    if temp_dir and temp_dir.exists():
        shutil.rmtree(temp_dir)
        logger.info(f"Cleaned up temporary directory {temp_dir}")
    
    if zip_path and zip_path.exists():
        zip_path.unlink()
        logger.info(f"Cleaned up zip file {zip_path}")


def main():
    parser = argparse.ArgumentParser(description="Update MSigDB gene set libraries")
    parser.add_argument(
        "--zip-file",
        type=Path,
        help="Path to MSigDB zip file"
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download MSigDB files directly (requires login)"
    )
    parser.add_argument(
        "--email",
        type=str,
        help="Email for MSigDB login (required with --download)"
    )
    parser.add_argument(
        "--password",
        type=str,
        help="Password for MSigDB login (required with --download)"
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip creating backup of existing libraries"
    )
    parser.add_argument(
        "--skip-background-update",
        action="store_true",
        help="Skip updating background gene lists"
    )
    parser.add_argument(
        "--skip-ncbi-update",
        action="store_true",
        help="Skip downloading NCBI gene information files"
    )
    
    args = parser.parse_args()
    
    if args.download:
        if not args.email or not args.password:
            logger.error("Email and password are required for direct download")
            return 1
    elif not args.zip_file:
        logger.error("Either --zip-file or --download must be specified")
        return 1
    
    try:
        # Download NCBI files (unless skipped)
        gene_info_path = None
        gene_history_path = None
        if not args.skip_ncbi_update:
            logger.info("Starting NCBI gene information update...")
            gene_info_gz_path, gene_history_gz_path = download_ncbi_files()
            
            if gene_info_gz_path and gene_history_gz_path:
                gene_info_path, gene_history_path = extract_gene_info_files(gene_info_gz_path, gene_history_gz_path)
            elif gene_info_gz_path:
                gene_info_path, _ = extract_gene_info_files(gene_info_gz_path, None)
            
            if gene_info_path:
                # Copy to project root for compatibility
                project_gene_info = PROJECT_ROOT / "Homo_sapiens.gene_info"
                shutil.copy2(gene_info_path, project_gene_info)
                logger.info(f"Copied gene_info to {project_gene_info}")
        else:
            logger.info("Skipping NCBI gene information update as requested")
            # Check if files exist in recent_release
            gene_info_path = RECENT_RELEASE_DIR / "Homo_sapiens.gene_info"
            gene_history_path = RECENT_RELEASE_DIR / "gene_history"
            
            # Check if compressed files exist and extract them
            gene_info_gz_path = RECENT_RELEASE_DIR / "Homo_sapiens.gene_info.gz"
            gene_history_gz_path = RECENT_RELEASE_DIR / "gene_history.gz"
            
            if not gene_info_path.exists() and gene_info_gz_path.exists():
                logger.info("Extracting existing gene_info.gz...")
                gene_info_path, _ = extract_gene_info_files(gene_info_gz_path, None)
            
            if not gene_history_path.exists() and gene_history_gz_path.exists():
                logger.info("Extracting existing gene_history.gz...")
                _, gene_history_path = extract_gene_info_files(None, gene_history_gz_path)
            
            if not gene_info_path.exists():
                logger.warning("Homo_sapiens.gene_info not found in recent_release")
                gene_info_path = None
            if not gene_history_path.exists():
                logger.warning("gene_history not found in recent_release")
                gene_history_path = None
        
        # Create backup unless disabled
        if not args.no_backup:
            create_backup()
        
        # Get the zip file
        if args.download:
            zip_path = download_msigdb(args.email, args.password)
            if not zip_path:
                return 1
        else:
            zip_path = args.zip_file
            if not zip_path.exists():
                logger.error(f"Zip file not found: {zip_path}")
                return 1
        
        # Extract the zip file
        temp_dir = extract_zip_file(zip_path)
        
        # Find .gmt files
        gmt_files = find_gmt_files(temp_dir)
        
        if not gmt_files:
            logger.error("No .gmt files found in the zip file")
            cleanup_temp_files(temp_dir, zip_path if args.download else None)
            return 1
        
        # Update libraries
        update_libraries(gmt_files)
        
        # Update alias.json
        update_alias_json()
        
        # Update background gene lists (unless skipped)
        if not args.skip_background_update:
            logger.info("Starting background gene list update for C2 and C5 gene sets (excluding c2.all and c5.all)...")
            
            # Collect all genes from C2 and C5 GMT files
            all_genes = collect_all_genes_from_gmt_files(gmt_files)
            
            # Validate and map genes
            if gene_info_path and gene_history_path:
                valid_genes, invalid_genes = validate_and_map_genes(all_genes, gene_info_path, gene_history_path)
            elif gene_info_path:
                valid_genes, invalid_genes = validate_and_map_genes(all_genes, gene_info_path, None)
            else:
                logger.warning("No gene_info available, skipping validation")
                valid_genes = all_genes
                invalid_genes = set()
            
            # Update background gene lists
            if valid_genes:
                update_background_gene_lists(valid_genes)
            else:
                logger.warning("No valid genes found, skipping background update")
        else:
            logger.info("Skipping background gene list update as requested")
        
        # Cleanup
        cleanup_temp_files(temp_dir, zip_path if args.download else None)
        
        logger.info("MSigDB update completed successfully!")
        return 0
        
    except Exception as e:
        logger.error(f"Error during update: {e}")
        return 1


if __name__ == "__main__":
    import sys
    exit(main())
