#!/bin/bash


set -euo pipefail

echo "🧬 Downloading HTLV-1 AF033817_1 for RNAseq mapping"

# Define file URLs
fasta_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/863/585/GCF_000863585.1_ViralProj15434/GCF_000863585.1_ViralProj15434_genomic.fna.gz"
gtf_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/863/585/GCF_000863585.1_ViralProj15434/GCF_000863585.1_ViralProj15434_genomic.gtf.gz"
id_url="AF033817_1"

# Define output directory
ref_dir="/home/htlvatl/coreapps/2025_refgenomes"

# Create output directory if it doesn't exist
mkdir -p "$ref_dir"

# Check if directory is writable
if [ ! -w "$ref_dir" ]; then
    echo "❌ Error: Directory '$ref_dir' is not writable"
    exit 1
fi

# Function to download files
download_file() {
    local url="$1"
    local output_path="$2"
    local description="$3"
    
    echo "📥 Downloading $description..."
    if wget -q --show-progress "$url" -O "$output_path"; then
        echo "✅ $description downloaded successfully"
        return 0
    else
        echo "❌ Error downloading $description"
        return 1
    fi
}


# Download FASTA file
download_file "$fasta_url" "$ref_dir/${id_url}.fa.gz" "HTLV-1 FASTA" || exit 1

# Download GTF file
download_file "$gtf_url" "$ref_dir/${id_url}.gtf.gz" "HTLV-1 GTF annotation" || exit 1


# Decompress files
echo "🔧 Decompressing files..."
cd "$ref_dir"

fasta_gz="${id_url}.fa.gz"
gtf_gz="${id_url}.gtf.gz"

if [ -f "$fasta_gz" ]; then
    echo "📦 Decompressing FASTA..."
    gunzip -f "$fasta_gz"
    echo "✅ FASTA decompressed: ${id_url}.fa"
else
    echo "❌ FASTA file not found: $fasta_gz"
    exit 1
fi

if [ -f "$gtf_gz" ]; then
    echo "📦 Decompressing GTF..."
    gunzip -f "$gtf_gz"
    echo "✅ GTF decompressed: ${id_url}.gtf"
else
    echo "❌ GTF file not found: $gtf_gz"
    exit 1
fi
