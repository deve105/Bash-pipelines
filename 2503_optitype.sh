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
	mkdir -p "${filename}"/optitype
else
	echo "New folder '$filename' created."
	mkdir -p "${filename}"/optitype #Creates the folder structure
fi

# Previously bwa-mem index should be configured.
if ! command -v samtools &>/dev/null; then
    printf "Error: samtools/OptiTypePipeline is either not installed or not in the PATH.\n"
    exit 1
fi

project_directory="/home/labatl/devprojects/Peru_IRID/optitype"
# Step 2 - Process each SRA identifier
while IFS= read -r sra; do
    echo "HLA typing for ${sra}"

    # Define paths
    bamfile="${filename}/bam/${sra}.bam"
    bamoutput="${filename}/optitype/${sra}"

    # Check if the BAM file exists
    if [ ! -f "$bamfile" ]; then
        echo "Error: No BAM file found for $sra at $bamfile. Skipping."
        continue
    fi

    # Step 3 - Extract reads mapped to Chromosome 6 (HLA)
    echo "Extracting reads mapped to Chromosome 6 (HLA) for ${sra}"
    samtools view -b "$bamfile" chr6 > "${bamoutput}_chr6.bam"

    # Convert BAM to FASTQ
    echo "Converting BAM to FASTQ"
    samtools sort -n "${bamoutput}_chr6.bam" -o "${bamoutput}_sorted_chr6.bam" 2>/dev/null
    bedtools bamtofastq -i "${bamoutput}_sorted_chr6.bam" -fq "${bamoutput}_sorted_chr6_1.fastq" 2>/dev/null -fq2 "${bamoutput}_sorted_chr6_2.fastq"

    # Clean up temporary files
    rm "${bamoutput}_chr6.bam" "${bamoutput}_sorted_chr6.bam"

    fastq_1="${sra}_sorted_chr6_1.fastq"
    fastq_2="${sra}_sorted_chr6_2.fastq"
    sudo docker run \
        -v "${project_directory}:/data/" \
        -t fred2/optitype \
        -d \
        -i "/data/${fastq_1}" "/data/${fastq_2}" \
        -o "/data/output/" \
        -p "${sra}_"

done < "$1"

echo "End of the script"