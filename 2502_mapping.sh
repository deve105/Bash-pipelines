#!/bin/bash

set -euo pipefail #Stric mode

echo "Script for SRA download and QC by Daniel Enriquez-Vera"
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

# Step 1 verify the folder exists
filename=$(basename "$1")
filename="${filename%.*}" # Removes the last extension

echo "Step 1. verifying that folder ${filename}/rawdata exists"
if [ -d "${filename}/rawdata" ] && [ "$(ls -A ${filename}/rawdata)" ]; then
	echo "Folder '${filename}/rawdata' detected"
else
	printf  "Error: '%s' a non-empty folder is needed.\n" "${filename}/rawdata"
    exit 1
fi

# Previously bwa-mem index should be configured.
if ! command -v bwa-mem2 &>/dev/null; then
    printf "Error: bwa-mem2 is either not installed or not in the PATH.\n"
    exit 1
fi

# Step 2 - Check if bwa-mem2 is installed
while IFS= read -r sra; do
    echo "Mapping ${sra}"

indexgenome="/home/labatl/devapps/2407_references/2405project_hg38/2405_hg38htlv.fa"


# Step 3 - Check if single-end or pair-end
if [ -f "${filename}/rawdata/${sra}_postqc_1.fq.gz" ] && [ -f "${filename}/rawdata/${sra}_postqc_2.fastq" ]; then
    # Paired-end reads
    fastq1="${filename}/rawdata/${sra}_postqc_1.fq.gz"
    fastq2="${filename}/rawdata/${sra}_postqc_2.fq.gz"
    echo "Paired-end reads detected for $sra"
    bwa-mem2 mem -t 20 "${indexgenome}" "${fastq1}" "${fastq2}" | \
    samtools view -o "${filename}/rawdata/${sra}.bam"

elif [ -f "${filename}/rawdata/${sra}_postqc_1.fq.gz" ]; then
    # Single-end reads
    fastq1="${filename}/rawdata/${sra}_postqc_1.fq.gz"
    echo "Single-end reads detected for $sra"
    bwa-mem2 mem -t 20 "${indexgenome}" "${fastq1}" | \
    samtools view -o "${filename}/rawdata/${sra}.bam"
else
    echo "Error: No valid FASTQ files found for $sra"
    exit 1
fi

done < "$1"

echo "End of the script"
