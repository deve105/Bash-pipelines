#!/bin/bash

set -euo pipefail

echo "Downloading GRCh38 reference genome for HLA analysis"
echo "Date: $(date)"

# Create reference directory
REF_DIR="/home/labatl/devapps/2407_references/GRCh38"
mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

# Download the best reference for HLA analysis
echo "Downloading GRCh38 with HLA decoy sequences..."
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Verify download
if [ -f "GRCh38_full_analysis_set_plus_decoy_hla.fa" ]; then
    echo "Download successful!"
    echo "File size: $(ls -lh GRCh38_full_analysis_set_plus_decoy_hla.fa | awk '{print $5}')"
    
    # Create BWA-MEM2 index
    echo "Creating BWA-MEM2 index..."
    echo "This will take 30-60 minutes..."
    
    # Use your BWA-MEM2 path
    BWA_MEM2="/home/htlvatl/Documents/devapps/bwa-mem2-2.3_x64-linux/bwa-mem2.avx2"
    
    if [ -f "${BWA_MEM2}" ]; then
        "${BWA_MEM2}" index GRCh38_full_analysis_set_plus_decoy_hla.fa
        echo "BWA-MEM2 index created successfully!"
    else
        echo "BWA-MEM2 not found at ${BWA_MEM2}"
        echo "Please update the path or install BWA-MEM2"
    fi
    
    # Create samtools index
    echo "Creating samtools faidx index..."
    samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
    
    echo "Reference genome setup complete!"
    echo "Location: ${PWD}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
else
    echo "Download failed!"
    exit 1
fi
