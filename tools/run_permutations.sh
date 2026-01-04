#!/bin/bash
#
# Run permutation distribution generation for iGEA analysis
# 
# This script runs the permutation distribution generator on Code Ocean.
# It automatically detects CPU cores and can be customized via environment variables.
#
# Usage:
#   ./run_permutations.sh
#   ./run_permutations.sh --n-permutations 500 --sizes 50 100 200
#
# Environment variables:
#   N_PERMUTATIONS: Number of permutations per size (default: 1000)
#   N_JOBS: Number of parallel jobs (default: auto-detect CPU cores)
#   RESUME: Resume from existing results (default: true, set to "false" to disable)
#   SIZES: Space-separated list of gene list sizes (default: all sizes 50-1000)
#
# Example for c6i.8xlarge (32 cores):
#   N_JOBS=32 ./run_permutations.sh
#

set -e  # Exit on error
set -u  # Exit on undefined variable

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Change to project root directory
cd "$PROJECT_ROOT"

# Default values
DEFAULT_N_PERMUTATIONS=1000
DEFAULT_N_JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo "4")
DEFAULT_RESUME="true"

# Get values from environment variables or use defaults
N_PERMUTATIONS=${N_PERMUTATIONS:-$DEFAULT_N_PERMUTATIONS}
N_JOBS=${N_JOBS:-$DEFAULT_N_JOBS}
RESUME=${RESUME:-$DEFAULT_RESUME}

# Build command arguments
ARGS=()

# Add n-permutations if specified
if [ -n "${N_PERMUTATIONS:-}" ]; then
    ARGS+=("--n-permutations" "$N_PERMUTATIONS")
fi

# Add n-jobs if specified
if [ -n "${N_JOBS:-}" ]; then
    ARGS+=("--n-jobs" "$N_JOBS")
fi

# Add resume flag
if [ "$RESUME" = "false" ] || [ "$RESUME" = "0" ]; then
    ARGS+=("--no-resume")
fi

# Add sizes if specified
if [ -n "${SIZES:-}" ]; then
    ARGS+=("--sizes" $SIZES)
fi

# Add any additional command-line arguments
ARGS+=("$@")

# Print configuration
echo "=========================================="
echo "iGEA Permutation Distribution Generator"
echo "=========================================="
echo "Project root: $PROJECT_ROOT"
echo "Script location: $SCRIPT_DIR"
echo "Number of permutations per size: $N_PERMUTATIONS"
echo "Number of parallel jobs: $N_JOBS"
echo "Resume mode: $RESUME"
if [ -n "${SIZES:-}" ]; then
    echo "Gene list sizes: $SIZES"
else
    echo "Gene list sizes: All (50-1000 in increments of 50)"
fi
echo "Output directory: results/permutation_results/"
echo "=========================================="
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null && ! command -v python &> /dev/null; then
    echo "Error: Python not found. Please ensure Python is installed."
    exit 1
fi

# Use python3 if available, otherwise python
PYTHON_CMD=$(command -v python3 2>/dev/null || command -v python 2>/dev/null)

# Check if the permutation script exists
SCRIPT_PATH="$SCRIPT_DIR/generate_permutation_distribution.py"
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: Permutation script not found at $SCRIPT_PATH"
    exit 1
fi

# Create results directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/results/permutation_results"

# Run the permutation script
echo "Starting permutation distribution generation..."
echo "Command: $PYTHON_CMD $SCRIPT_PATH ${ARGS[*]}"
echo ""

# Run with error handling
if $PYTHON_CMD "$SCRIPT_PATH" "${ARGS[@]}"; then
    echo ""
    echo "=========================================="
    echo "Permutation distribution generation completed successfully!"
    echo "=========================================="
    echo "Results saved to: results/permutation_results/"
    echo "Summary file: results/permutation_results/summary.json"
    echo "Config file: results/permutation_results/config.json"
    echo ""
    exit 0
else
    EXIT_CODE=$?
    echo ""
    echo "=========================================="
    echo "Error: Permutation distribution generation failed with exit code $EXIT_CODE"
    echo "=========================================="
    echo "Check the log file: permutation_distribution.log"
    echo ""
    exit $EXIT_CODE
fi

