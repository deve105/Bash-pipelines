#!/bin/bash

set -euo pipefail #Stric mode

echo "subsetting for optitype by Daniel Enriquez-Vera"
echo "Version V1.2 2025-02-20"

#-If an argument is needed
if [ $# -eq 0 ]; then
	printf "Usage: %s project.txt is needed.\n" "$0"
	exit 1
fi

#-If an existing file is needed
if [ ! -f "$1" ]; then
	printf  "Error: '%s' A valid file is needed.\n" "$1"
	exit 1
fi

#-If a non-empty file is needed 
if [ ! -s "$1" ]; then
	printf  "Error: '%s' a non-empty file is needed.\n" "$1"
	exit 1
fi

# Step 1 - create a subfolder with the project ID
filename=$(basename "$1")
filename="${filename%.*}" # Removes the last extension

echo "Step 1. Creating e a folder for the project '$filename'"
if [ -d "$filename" ]; then
	echo "Folder '$filename' already exists. Be careful"
	mkdir -p "${filename}/HLA-LA/output"
else
	echo "New folder '$filename' created."
	mkdir -p "${filename}/HLA-LA" #Creates the folder structure
fi

# Previously bwa-mem index should be configured.
if ! command -v HLA-LA &>/dev/null; then
    printf "Error: HLA-LA is either not installed or not in the PATH.\n"
    exit 1
fi

project_directory="/home/labatl/devprojects/Peru_IRID/HLA-LA"
# Step 2 - Process each SRA identifier
while IFS= read -r sra; do
    echo "HLA typing for ${sra}"

    # Define paths
    bamfile="${filename}/bam/${sra}.bam"
    #bamoutput="${filename}/optitype/${sra}"

    # Check if the BAM file exists
    if [ ! -f "$bamfile" ]; then
        echo "Error: No BAM file found for $sra at $bamfile. Skipping."
        continue
    fi

    # Step 3 - Extract reads mapped to Chromosome 6 (HLA)
    echo "Extracting reads mapped to Chromosome 6 (HLA) for ${sra}"
    HLA-LA.pl --BAM "${bamfile}" \
        --graph PRG_MHC_GRCh38_withIMGT \
        --sampleID "${sra}" \
        --maxThreads 20 \
        --workingDir "${project_directory}" 
done < "$1"
echo "End of the script"