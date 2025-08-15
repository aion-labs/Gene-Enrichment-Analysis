#!/bin/bash

# Gene Enrichment Analysis - Linux Installation Script
# This script automates the installation process on a clean Linux machine

set -e  # Exit on any error

echo "🧬 Gene Enrichment Analysis - Linux Installation"
echo "================================================"

# Check if running as root
if [[ $EUID -eq 0 ]]; then
   echo "❌ This script should not be run as root. Please run as a regular user."
   exit 1
fi

# Detect Linux distribution
if [ -f /etc/os-release ]; then
    . /etc/os-release
    OS=$NAME
    VER=$VERSION_ID
else
    echo "❌ Cannot detect Linux distribution"
    exit 1
fi

echo "📋 Detected OS: $OS $VER"

# Function to install packages based on distribution
install_system_packages() {
    echo "📦 Installing system packages..."
    
    if [[ "$OS" == *"Ubuntu"* ]] || [[ "$OS" == *"Debian"* ]]; then
        sudo apt update
        sudo apt install -y python3 python3-pip python3-venv python3-dev \
            build-essential gfortran libopenblas-dev liblapack-dev \
            graphviz graphviz-dev libffi-dev libssl-dev zlib1g-dev \
            libbz2-dev libreadline-dev libsqlite3-dev git wget curl
    elif [[ "$OS" == *"CentOS"* ]] || [[ "$OS" == *"Red Hat"* ]]; then
        sudo yum update -y
        sudo yum groupinstall -y "Development Tools"
        sudo yum install -y python3 python3-pip python3-devel \
            gcc-gfortran openblas-devel lapack-devel graphviz graphviz-devel \
            libffi-devel openssl-devel zlib-devel bzip2-devel readline-devel \
            sqlite-devel git wget curl
    elif [[ "$OS" == *"Fedora"* ]]; then
        sudo dnf update -y
        sudo dnf groupinstall -y "Development Tools"
        sudo dnf install -y python3 python3-pip python3-devel \
            gcc-gfortran openblas-devel lapack-devel graphviz graphviz-devel \
            libffi-devel openssl-devel zlib-devel bzip2-devel readline-devel \
            sqlite-devel git wget curl
    else
        echo "❌ Unsupported Linux distribution: $OS"
        echo "Please install the required packages manually and run this script again."
        exit 1
    fi
}

# Function to create virtual environment
setup_virtual_environment() {
    echo "🐍 Setting up Python virtual environment..."
    
    if [ ! -d ".venv" ]; then
        python3 -m venv .venv
    fi
    
    source .venv/bin/activate
    pip install --upgrade pip
}

# Function to install Python dependencies
install_python_dependencies() {
    echo "📚 Installing Python dependencies..."
    
    if [ -f "requirements.txt" ]; then
        pip install -r requirements.txt
    else
        echo "❌ requirements.txt not found. Installing core dependencies manually..."
        pip install pandas==2.2.3 scipy==1.15.2 statsmodels==0.14.4 \
            plotly==6.0.1 pydot==4.0.1 networkx==3.5 streamlit==1.43.2 \
            typer==0.9.0 Pillow==11.1.0 graphviz==0.21 numpy matplotlib seaborn
    fi
}

# Function to download data files
download_data_files() {
    echo "📥 Downloading required data files..."
    
    # Create directories
    mkdir -p data/backgrounds data/libraries data/recent_release
    
    # Download NCBI gene information
    if [ ! -f "Homo_sapiens.gene_info" ]; then
        echo "📥 Downloading NCBI gene information..."
        wget -O Homo_sapiens.gene_info.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
        gunzip Homo_sapiens.gene_info.gz
    else
        echo "✅ NCBI gene information already exists"
    fi
    
    # Download gene history
    if [ ! -f "data/recent_release/gene_history" ]; then
        echo "📥 Downloading gene history..."
        wget -O data/recent_release/gene_history.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_HISTORY.gz
        gunzip data/recent_release/gene_history.gz
    else
        echo "✅ Gene history already exists"
    fi
    
    # Download background gene set (if available)
    if [ ! -f "data/backgrounds/all_genes.txt" ]; then
        echo "📥 Downloading background gene set..."
        wget -O data/backgrounds/all_genes.txt https://raw.githubusercontent.com/aion-labs/Gene-Enrichment-Analysis/main/data/backgrounds/all_genes.txt 2>/dev/null || echo "⚠️  Could not download background gene set. Please download manually."
    else
        echo "✅ Background gene set already exists"
    fi
}

# Function to test installation
test_installation() {
    echo "🧪 Testing installation..."
    
    # Test Python imports
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
    
    # Test app modules
    python -c "
import sys
sys.path.append('code')
try:
    from code.gene_set import GeneSet
    from code.gene_set_library import GeneSetLibrary
    from code.background_gene_set import BackgroundGeneSet
    print('✅ Core app modules imported successfully')
except ImportError as e:
    print(f'⚠️  Some modules could not be imported: {e}')
"
}

# Function to create launcher script
create_launcher() {
    echo "🚀 Creating launcher script..."
    
    cat > launch_gene_enrichment.sh << 'EOF'
#!/bin/bash

# Gene Enrichment Analysis Launcher Script
# This script launches the Streamlit app for gene enrichment analysis

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Change to the project directory
cd "$SCRIPT_DIR"

# Check if virtual environment exists
if [ ! -d ".venv" ]; then
    echo "❌ Virtual environment not found. Please run the installation script first."
    exit 1
fi

# Activate virtual environment
echo "🔧 Activating virtual environment..."
source .venv/bin/activate

# Check if required packages are installed
if ! python -c "import streamlit" 2>/dev/null; then
    echo "❌ Streamlit not found. Installing requirements..."
    pip install -r requirements.txt 2>/dev/null || pip install streamlit plotly pandas scipy statsmodels
fi

# Launch the Streamlit app
echo "🚀 Launching Gene Enrichment Analysis..."
echo "📊 App will open in your browser at: http://localhost:8501"
echo "🛑 Press Ctrl+C to stop the app"
echo ""

streamlit run code/streamlit_app.py
EOF

    chmod +x launch_gene_enrichment.sh
}

# Function to display next steps
display_next_steps() {
    echo ""
    echo "🎉 Installation completed successfully!"
    echo "======================================"
    echo ""
    echo "📋 Next steps:"
    echo "1. Download MSigDB libraries from https://www.gsea-msigdb.org/gsea/msigdb/"
    echo "2. Place the .gmt files in the data/libraries/ directory"
    echo "3. Run the library validator: python update_msigdb_libraries.py data/libraries/"
    echo "4. Launch the app: ./launch_gene_enrichment.sh"
    echo ""
    echo "🌐 The app will be available at: http://localhost:8501"
    echo ""
    echo "📚 For more information, see LINUX_INSTALLATION.md"
}

# Main installation process
main() {
    echo "🚀 Starting installation process..."
    
    # Install system packages
    install_system_packages
    
    # Setup virtual environment
    setup_virtual_environment
    
    # Install Python dependencies
    install_python_dependencies
    
    # Download data files
    download_data_files
    
    # Test installation
    test_installation
    
    # Create launcher script
    create_launcher
    
    # Display next steps
    display_next_steps
}

# Run main function
main

