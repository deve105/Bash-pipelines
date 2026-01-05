#!/bin/bash
# filepath: /home/htlvatl/Documents/GitHub/Bash-pipelines/2510_RNAseq/2601_041_sc_download_index.sh

set -euo pipefail

echo "🧬 Downloading Human Index Genome GRCh38 v49 for RNAseq mapping"

# Define file URLs
fasta_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"
gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz"

# Define output directory
ref_dir="/home/htlvatl/coreapps/2025_refgenomes"

# Create output directory if it doesn't exist
mkdir -p "$ref_dir"

# Check if directory is writable
if [ ! -w "$ref_dir" ]; then
    echo "❌ Error: Directory '$ref_dir' is not writable"
    exit 1
fi

# Download FASTA file
echo "📥 Downloading FASTA file..."
if wget -P "$ref_dir" "$fasta_url"; then
    echo "✅ FASTA file downloaded successfully"
else
    echo "❌ Error downloading FASTA file"
    exit 1
fi

# Download GTF file
echo "📥 Downloading GTF annotation file..."
if wget -P "$ref_dir" "$gtf_url"; then
    echo "✅ GTF file downloaded successfully"
else
    echo "❌ Error downloading GTF file"
    exit 1
fi

# Decompress files
echo "🔧 Decompressing files..."
cd "$ref_dir"

fasta_gz="GRCh38.primary_assembly.genome.fa.gz"
gtf_gz="gencode.v49.primary_assembly.annotation.gtf.gz"

if [ -f "$fasta_gz" ]; then
    echo "📦 Decompressing FASTA..."
    gunzip -f "$fasta_gz"
    echo "✅ FASTA decompressed"
else
    echo "❌ FASTA file not found"
    exit 1
fi

if [ -f "$gtf_gz" ]; then
    echo "📦 Decompressing GTF..."
    gunzip -f "$gtf_gz"
    echo "✅ GTF decompressed"
else
    echo "❌ GTF file not found"
    exit 1
fi
