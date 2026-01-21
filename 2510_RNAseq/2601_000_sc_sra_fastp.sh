#!/bin/bash

set -euo pipefail #Stric mode

echo "Script for SRA download, fastp QC, STAR mapping, and multiQC" 
echo "by Daniel Enriquez Version V1.2 2026-01-07"


disc_tmp="/home/labatl/devprojects"

filename=$(basename "$1")
filename="${filename%.*}" # Removes the last extension

filename_disc="${disc_tmp}/${filename}"
platform="ILLUMINA"
numthreads=24




## ============================================================================
## PRE-REQUISITES CHECKS - File arguments
## ============================================================================
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


## ============================================================================
## PRE-REQUISITES CHECKS - Folder arguments
## ============================================================================



# Step 2 create temp project folder
echo "Step 1. Create a Temporal folder for the project in HDD'$filename'"
if [ -d "$filename_disc" ]; then
	echo "Folder '$filename_disc' already exists. Be careful"
	mkdir -p "${filename_disc}"/{rawdata,bam,logs,reports}
else
	echo "New folder '$filename_disc' created."
	mkdir -p "${filename_disc}"/{rawdata,bam,logs,reports} #Creates the folder structure
fi



# ============================================================================
# STEP 1 - CHECK DEPENDENCIES
# ============================================================================

echo "🔍 Checking dependencies..."

for tool in prefetch fasterq-dump fastp samtools multiqc; do
    if ! command -v "$tool" &>/dev/null; then
        printf "❌ Error: %s is not installed or not in PATH\n" "$tool"
        exit 1
    fi
done

echo "✅ All dependencies found"



# ============================================================================
# STEP 2 - PROCESS SRA IDS
# ============================================================================


# Step 3 - Download SRA data
while IFS= read -r sra; do
    cd "${disc_tmp}"
    
    echo "Downloading ${sra}"

    # Step 3.1 - Prefetch the SRA data
    if ! prefetch "${sra}"; then
        echo "Error: prefetch failed for ${sra}"
        exit 1
    fi

    # Step 3.2 - Convert SRA to FASTQ using fasterq-dump
    if ! fasterq-dump --mem 30G \
                      --threads "${numthreads}" \
                      --skip-technical \
                      --outdir "${filename_disc}/rawdata/" \
                      --temp "${filename_disc}" \
                      --split-3 \
                      "${sra}" \
                      2>&1 | tee -a "${filename_disc}/logs/${sra}_fasterq_dump.log"; then
        echo "Error: fasterq-dump failed for ${sra}"
        
        exit 1
    fi
    rm -rf "${disc_tmp}/${sra}"
		
	# Check if the files are paired-end or single-end
    if [ -f "${filename_disc}/rawdata/${sra}_1.fastq" ] && [ -f "${filename_disc}/rawdata/${sra}_2.fastq" ]; then
        # Paired-end reads
        echo "Paired-end reads detected for ${sra}"
        gzip -f "${filename_disc}/rawdata/${sra}_1.fastq" "${filename_disc}/rawdata/${sra}_2.fastq"
    elif [ -f "${filename_disc}/rawdata/${sra}.fastq" ]; then
        # Single-end reads
        echo "Single-end reads detected for $sra"
        mv "${filename_disc}/rawdata/${sra}.fastq" "${filename_disc}/rawdata/${sra}_1.fastq"
        gzip -f "${filename_disc}/rawdata/${sra}_1.fastq"
    else
        echo "Error: No valid FASTQ files found for $sra"
        exit 1
    fi
    echo "${sra} was downloaded and processed :)"
    if [ -f "${filename_disc}/rawdata/${sra}_1.fastq.gz" ] && [ -f "${filename_disc}/rawdata/${sra}_2.fastq.gz" ]; then
        echo "Performing fastp quality control for paired-end reads: ${sra}"
    # Paired-end mode
        fastp -i "${filename_disc}/rawdata/${sra}_1.fastq.gz" \
		    -I "${filename_disc}/rawdata/${sra}_2.fastq.gz" \
            -o "${filename_disc}/rawdata/${sra}_postqc_1.fq.gz" \
            -O "${filename_disc}/rawdata/${sra}_postqc_2.fq.gz" \
			-j "${filename_disc}/reports/${sra}_fastp.json" \
            -h "${filename_disc}/reports/${sra}_fastp.html" \
            --detect_adapter_for_pe \
            --length_required 15 \
            --correction \
            -e 15 \
            -q 15 \
            -w 12 \
            -u 40 \
            -n 5 \
            2>&1 | tee -a "${filename_disc}/logs/${sra}_fastp.log"

    elif [ -f "${filename_disc}/rawdata/${sra}_1.fastq.gz" ]; then
        echo "Performing fastp quality control for single-end reads: ${sra}"
        # Single-end mode
        fastp -i "${filename_disc}/rawdata/${sra}_1.fastq.gz" \
            -o "${filename_disc}/rawdata/${sra}_postqc_1.fq.gz" \
			-j "${filename_disc}/reports/${sra}_fastp.json" \
            -h "${filename_disc}/reports/${sra}_fastp.html" \
            --length_required 15 \
            --correction \
            -e 15 \
            -q 15 \
            -w 12 \
            -u 40 \
            -n 5 \
            2>&1 | tee -a "${filename_disc}/logs/${sra}_fastp.log"
    else
        echo "Error: No valid FASTQ files found for ${sra}"
        exit 1
    fi
    echo "Fastp QC completed for ${sra}"

done < "$1"