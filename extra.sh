#!/bin/bash

set -euo pipefail #Stric mode

echo "Script for SRA download by Daniel Enriquez-Vera"
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

# Step 1 - create a folder with the project ID
filename=$(basename "$1")
filename="${filename%.*}" # Removes the last extension

echo "Step 1. Create a folder for the project '$filename'"
if [ -d "$filename" ]; then
	echo "Folder '$filename' already exists. Be careful"
	mkdir -p "${filename}"/{rawdata,bam,reports}
else
	echo "New folder '$filename' created."
	mkdir -p "${filename}"/{rawdata,bam,reports} #Creates the folder structure
fi

# Step 2 - Check if sratoolkit is installed
# Previously vdb-config -i should be configured.
if ! command -v prefetch &>/dev/null || ! command -v fasterq-dump &>/dev/null; then
    printf "Error: sratools is either not installed or not in the PATH.\n"
    exit 1
fi

if ! command -v multiqc &>/dev/null; then
    printf "Error: MultiQC is not installed or not in the PATH.\n"
    exit 1
fi


echo "Running fastp on the rawdata folder"
while IFS= read -r sra; do
    if [ -f "${filename}/rawdata/${sra}_1.fq.gz" ] && [ -f "${filename}/rawdata/${sra}_2.1.fq.gz" ]; then
    # Paired-end mode
        fastp -i "${filename}/rawdata/${sra}_1.fq.gz" \
		    -I "${filename}/rawdata/${sra}_2.1.fq.gz" \
            -3 \
            -o "${filename}/rawdata/${sra}_postqc_1.fq.gz" \
            -O "${filename}/rawdata/${sra}_postqc_2.fq.gz" \
			-j "${filename}/reports/${sra}_fastp.json" \
            -h "${filename}/reports/${sra}_fastp.html" \
            -q 20 \ #Sets the quality threshold for filtering reads. Reads with average quality below this 
            -u 40 \ #Sets the percentage of bases allowed to be below the quality threshold
            -n 5 \ #Discards reads containing more than this number of N bases
		rm "${filename}/rawdata/${sra}_1.fq.gz" "${filename}/rawdata/${sra}_2.fq.gz"	
    elif [ -f "${filename}/rawdata/${sra}_1.fq.gz" ]; then
        # Single-end mode
        fastp -i "${filename}/rawdata/${sra}_1.fq.gz"  \
		    -o "${filename}/rawdata/${sra}_postqc_1.fq.gz" \
            -j "${filename}/reports/${sra}_fastp.json" \
			-h "${filename}/reports/${sra}_fastp.html" \
            -q 20 \
            -u 40 \
            -n 5 \ 
		rm "${filename}/rawdata/${sra}_1.fq.gz" 
    else
        echo "Error: No valid FASTQ files found for $sra"
        exit 1
    fi
done < "$1"

       HLA-LA.pl --BAM /home/labatl/devprojects/Peru_IRID/bam/SS02_sorted.bam \
        --graph PRG_MHC_GRCh38_withIMGT \
        --sampleID SS02 \
        --maxThreads 30 \
        --workingDir /home/labatl/devprojects/Peru_IRID/HLA-LA \
        --samtools_T /home/labatl/devapps/2407_references/2405project_hg38/2405_hg38htlv.fa  