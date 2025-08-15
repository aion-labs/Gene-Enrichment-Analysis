# Linux Installation Guide

This guide provides step-by-step instructions for installing and running the Gene Enrichment Analysis app on a clean Linux machine.

## System Requirements

- **OS**: Ubuntu 20.04+, CentOS 8+, or similar Linux distribution
- **Python**: 3.12.4 or higher
- **RAM**: Minimum 4GB, recommended 8GB+
- **Storage**: At least 2GB free space
- **Internet**: Required for downloading dependencies and data

## Installation Steps

### 1. System Package Installation

First, update your system and install required system packages:

```bash
# Update system packages
sudo apt update && sudo apt upgrade -y

# Install Python and development tools
sudo apt install -y python3 python3-pip python3-venv python3-dev

# Install system dependencies for scientific computing
sudo apt install -y build-essential gfortran libopenblas-dev liblapack-dev

# Install Graphviz (required for network visualization)
sudo apt install -y graphviz graphviz-dev

# Install additional system libraries
sudo apt install -y libffi-dev libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev

# Install Git (for cloning the repository)
sudo apt install -y git

# Install wget and curl (for downloading data)
sudo apt install -y wget curl
```

### 2. Clone the Repository

```bash
# Clone the repository
git clone https://github.com/aion-labs/Gene-Enrichment-Analysis.git
cd Gene-Enrichment-Analysis
```

### 3. Create Virtual Environment

```bash
# Create a virtual environment
python3 -m venv .venv

# Activate the virtual environment
source .venv/bin/activate
```

### 4. Install Python Dependencies

```bash
# Upgrade pip
pip install --upgrade pip

# Install core dependencies
pip install pandas==2.2.3
pip install scipy==1.15.2
pip install statsmodels==0.14.4
pip install plotly==6.0.1
pip install pydot==4.0.1
pip install networkx==3.5
pip install streamlit==1.43.2
pip install typer==0.9.0
pip install Pillow==11.1.0
pip install graphviz==0.21

# Install additional dependencies
pip install numpy
pip install matplotlib
pip install seaborn
```

### 5. Download Required Data Files

The app requires several data files to function properly:

```bash
# Create data directories
mkdir -p data/backgrounds data/libraries data/recent_release

# Download NCBI gene information (required for gene validation)
wget -O Homo_sapiens.gene_info.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
gunzip Homo_sapiens.gene_info.gz

# Download gene history (for gene symbol mapping)
wget -O data/recent_release/gene_history.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_HISTORY.gz
gunzip data/recent_release/gene_history.gz

# Download background gene set (all human genes)
wget -O data/backgrounds/all_genes.txt https://raw.githubusercontent.com/aion-labs/Gene-Enrichment-Analysis/main/data/backgrounds/all_genes.txt
```

### 6. Download MSigDB Libraries

You'll need to download the MSigDB gene set libraries. These are large files, so ensure you have sufficient bandwidth and storage:

```bash
# Create libraries directory
mkdir -p data/libraries

# Download MSigDB libraries (you'll need to get these from MSigDB website)
# For demonstration, here are the required files:
# - c1.all.v2025.1.Hs.symbols.gmt
# - c2.all.v2025.1.Hs.symbols.gmt
# - c2.cp.reactome.v2025.1.Hs.symbols.gmt
# - c5.go.bp.v2025.1.Hs.symbols.gmt
# - h.all.v2025.1.Hs.symbols.gmt
# - etc.

# Example (you'll need to replace with actual URLs):
# wget -O data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt [URL]
# wget -O data/libraries/c5.go.bp.v2025.1.Hs.symbols.gmt [URL]
```

**Note**: MSigDB requires registration. Visit [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) to download the required `.gmt` files.

### 7. Validate Libraries (Optional but Recommended)

```bash
# Run the library validator to ensure consistency
python update_msigdb_libraries.py data/libraries/
```

### 8. Test the Installation

```bash
# Test that all modules can be imported
python -c "
import pandas
import scipy
import statsmodels
import plotly
import streamlit
import networkx
import graphviz
print('✅ All core dependencies imported successfully')
"

# Test the app structure
python -c "
import sys
sys.path.append('code')
from code.gene_set import GeneSet
from code.gene_set_library import GeneSetLibrary
from code.background_gene_set import BackgroundGeneSet
print('✅ Core app modules imported successfully')
"
```

## Running the Application

### Method 1: Direct Streamlit Command

```bash
# Make sure virtual environment is activated
source .venv/bin/activate

# Run the Streamlit app
streamlit run code/streamlit_app.py
```

### Method 2: Using the Launcher Script

```bash
# Make the launcher script executable
chmod +x launch_gene_enrichment.sh

# Run the app
./launch_gene_enrichment.sh
```

### Method 3: Create a System Alias

```bash
# Add to your ~/.bashrc or ~/.zshrc
echo 'alias gene-enrichment="cd /path/to/Gene-Enrichment-Analysis && source .venv/bin/activate && streamlit run code/streamlit_app.py"' >> ~/.bashrc

# Reload shell configuration
source ~/.bashrc

# Use the alias
gene-enrichment
```

## Accessing the Application

Once running, the application will be available at:
- **Local**: http://localhost:8501
- **Network**: http://[your-ip]:8501

## Troubleshooting

### Common Issues

**1. Import Errors**
```bash
# If you get import errors, make sure you're in the virtual environment
source .venv/bin/activate
```

**2. Graphviz Issues**
```bash
# If graphviz fails, install system packages
sudo apt install -y graphviz graphviz-dev
```

**3. Memory Issues**
```bash
# If you encounter memory issues, increase swap space
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

**4. Permission Issues**
```bash
# Fix file permissions
chmod -R 755 data/
chmod -R 644 data/libraries/*.gmt
```

### Performance Optimization

**For Large Datasets:**
```bash
# Install optimized BLAS libraries
sudo apt install -y libopenblas-dev liblapack-dev

# Set environment variables for better performance
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
```

## Development Setup (Optional)

If you want to contribute to the project:

```bash
# Install development dependencies
pip install pytest==7.4.4
pip install isort==6.0.1
pip install black==25.1.0
pip install flake8==7.3.0

# Run tests
pytest

# Format code
black code/
isort code/
```

## System Requirements Summary

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| Python | 3.12.4 | 3.12.4+ |
| RAM | 4GB | 8GB+ |
| Storage | 2GB | 5GB+ |
| CPU | 2 cores | 4+ cores |

## Support

If you encounter issues:
1. Check the troubleshooting section above
2. Verify all dependencies are installed correctly
3. Ensure data files are in the correct locations
4. Check the application logs for error messages

For additional support, please refer to the main README.md file or create an issue on the GitHub repository.

