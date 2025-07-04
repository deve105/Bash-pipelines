#!/bin/bash

set -euo pipefail

echo "Updating BWA-MEM2 from source"
echo "Date: $(date)"

# Installation directory
INSTALL_DIR="/home/labatl/devapps/bwa-mem2"
BWA_MEM2_DIR="${INSTALL_DIR}/bwa-mem2"

# Check if existing installation exists
if [ -d "${BWA_MEM2_DIR}" ]; then
    echo "Found existing BWA-MEM2 installation at ${BWA_MEM2_DIR}"
    echo "Updating existing installation..."
    cd "${BWA_MEM2_DIR}"
else
    echo "No existing installation found. Creating new installation..."
    mkdir -p "${INSTALL_DIR}"
    cd "${INSTALL_DIR}"
fi

# Check if we have necessary build tools
echo "Checking for required build tools..."
if ! command -v gcc &> /dev/null; then
    echo "Error: gcc not found. Install with: sudo apt install build-essential"
    exit 1
fi

if ! command -v make &> /dev/null; then
    echo "Error: make not found. Install with: sudo apt install build-essential"
    exit 1
fi

if ! command -v git &> /dev/null; then
    echo "Error: git not found. Install with: sudo apt install git"
    exit 1
fi

# Handle updating or cloning
if [ -d "${BWA_MEM2_DIR}" ]; then
    echo "Updating existing repository..."
    cd "${BWA_MEM2_DIR}"
    
    # Clean any previous build artifacts
    echo "Cleaning previous build..."
    make clean 2>/dev/null || true
    
    # Fetch latest changes
    echo "Fetching latest changes..."
    git fetch origin
    
    # Check current branch and update
    CURRENT_BRANCH=$(git branch --show-current)
    echo "Current branch: ${CURRENT_BRANCH}"
    
    # Pull latest changes
    git pull origin "${CURRENT_BRANCH}"
    
    # Update submodules
    echo "Updating submodules..."
    git submodule update --init --recursive
else
    echo "Cloning BWA-MEM2 repository..."
    cd "${INSTALL_DIR}"
    git clone --recursive https://github.com/bwa-mem2/bwa-mem2
    cd bwa-mem2
    
    # Alternative method if recursive clone fails
    if [ ! -d "ext" ] || [ -z "$(ls -A ext 2>/dev/null)" ]; then
        echo "Initializing submodules..."
        git submodule init
        git submodule update
    fi
fi

# Check available memory before compilation
echo "Checking available memory..."
free -h

# Compile with optimizations and error handling
echo "Compiling BWA-MEM2..."
echo "This may take several minutes..."

# Backup existing binary if it exists
if [ -f "bwa-mem2" ]; then
    echo "Backing up existing binary..."
    cp bwa-mem2 bwa-mem2.backup.$(date +%Y%m%d_%H%M%S)
fi

# Use fewer parallel jobs to avoid memory issues
echo "Starting compilation (using $(nproc --ignore=2) cores)..."
if ! make -j$(nproc --ignore=2); then
    echo "Parallel compilation failed, trying single-threaded..."
    make clean
    make -j1
fi

# Check if compilation was successful
if [ -f "bwa-mem2" ]; then
    echo "Compilation successful!"
    echo "Testing the binary..."
    ./bwa-mem2
    
    # Update symlink in a directory that's likely in PATH
    if [ -L "/usr/local/bin/bwa-mem2" ] || [ -f "/usr/local/bin/bwa-mem2" ]; then
        echo "Updating existing symlink in /usr/local/bin..."
        sudo rm -f /usr/local/bin/bwa-mem2
    else
        echo "Creating new symlink in /usr/local/bin..."
    fi
    
    sudo ln -sf "${PWD}/bwa-mem2" /usr/local/bin/bwa-mem2
    
    echo "BWA-MEM2 updated successfully!"
    echo "Location: ${PWD}/bwa-mem2"
    echo "Symlink: /usr/local/bin/bwa-mem2"
    
    # Show git information
    echo ""
    echo "Git information:"
    echo "Commit: $(git rev-parse --short HEAD)"
    echo "Date: $(git log -1 --format=%cd --date=short)"
    echo "Branch: $(git branch --show-current)"
else
    echo "Compilation failed!"
    if [ -f "bwa-mem2.backup"* ]; then
        echo "Restoring backup..."
        cp bwa-mem2.backup.* bwa-mem2 2>/dev/null || true
    fi
    exit 1
fi

# Test the installation
echo ""
echo "=== Installation Summary ==="
echo "Testing installation..."
which bwa-mem2
echo ""
echo "Version information:"
bwa-mem2 version 2>/dev/null || echo "Version command not available"
echo ""
echo "Binary location: $(which bwa-mem2)"
echo "Binary size: $(ls -lh $(which bwa-mem2) | awk '{print $5}')"
echo "Installation complete!"
