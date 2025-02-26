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

# Step 3 - Download SRA data
while IFS= read -r sra; do
    echo "Downloading ${sra}"

    # Step 3.1 - Prefetch the SRA data
    if ! prefetch "${sra}"; then
        echo "Error: prefetch failed for ${sra}"
        exit 1
    fi

    # Step 3.2 - Convert SRA to FASTQ using fasterq-dump
    if ! fasterq-dump --mem 6G \
                      --threads 6 \
                      --skip-technical \
                      --outdir "${filename}/rawdata/" \
                      --temp "${filename}" \
                      --split-3 \
                      "${sra}"; then
        echo "Error: fasterq-dump failed for ${sra}"
        exit 1
    fi
		
	# Check if the files are paired-end or single-end
    if [ -f "${filename}/rawdata/${sra}_1.fastq" ] && [ -f "${filename}/rawdata/${sra}_2.fastq" ]; then
        # Paired-end reads
        echo "Paired-end reads detected for $sra"
        mv "${filename}/rawdata/${sra}_1.fastq" "${filename}/rawdata/${sra}_1.fq"
        mv "${filename}/rawdata/${sra}_2.fastq" "${filename}/rawdata/${sra}_2.fq"
        gzip "${filename}/rawdata/${sra}_1.fq" "${filename}/rawdata/${sra}_2.fq"
    elif [ -f "${filename}/rawdata/${sra}.fastq" ]; then
        # Single-end reads
        echo "Single-end reads detected for $sra"
        mv "${filename}/rawdata/${sra}.fastq" "${filename}/rawdata/${sra}_1.fq"
        gzip "${filename}/rawdata/${sra}_1.fq"
        # Create an empty _2.fq.gz file for consistency
        touch "${filename}/rawdata/${sra}_2.fq.gz"
    else
        echo "Error: No valid FASTQ files found for $sra"
        exit 1
    fi
    echo "${sra} was downloaded and processed :)"
done < "$1"

# Step 4 - Run MultiQC on the rawdata folder
echo "Running MultiQC on the rawdata folder"
multiqc "${filename}/rawdata/" -o "${filename}/reports/" -n "${filename}_multiqc_report_raw"

echo "MultiQC report generated in ${filename}/reports/${filename}_multiqc_report_raw.html"

# Step 5 - Run fastp on the rawdata folder
echo "Running fastp on the rawdata folder"
while IFS= read -r sra; do
    fastq1="${filename}/rawdata/${sra}_1.fq.gz"
    fastq2="${filename}/rawdata/${sra}_2.fq.gz"

    if [ -f "$fastq1" ] && [ -f "$fastq2" ]; then
        # Paired-end mode
        fastp -i "$fastq1" \
              -I "$fastq2" \
              -3 \
              -o "${filename}/rawdata/${sra}_postqc_1.fq.gz" \
              -O "${filename}/rawdata/${sra}_postqc_2.fq.gz" \
              -h "${filename}/reports/${sra}_fastp.html" \
              -q 20 \ #Sets the quality threshold for filtering reads. Reads with average quality below this 
              -u 40 \ #Sets the percentage of bases allowed to be below the quality threshold
              -n 5 \ #Discards reads containing more than this number of N bases
              -l 15 \ #Discards reads shorter than this length after trimming
              -w 20 \
              -t 0 \ #Trims a fixed number of bases from the 5-end of each read.
              -T 0 #Trims a fixed number of bases from the 3-end of each read.
    elif [ -f "$fastq1" ]; then
        # Single-end mode
        fastp -i "$fastq1" \
              -o "${filename}/rawdata/${sra}_postqc_1.fq.gz" \
              -h "${filename}/reports/${sra}_fastp.html" \
              -q 20 \
              -u 40 \
              -n 5 \ 
              -l 15 \ 
              -w 20 \
              -t 0 \
              -T 0
    else
        echo "Error: No valid FASTQ files found for $sra"
        exit 1
    fi

    # Remove the original files after processing
    rm "$fastq1" "$fastq2"
done < "$1"

# Step 6 - Re-Run MultiQC on the rawdata folder
echo "Running MultiQC on the rawdata folder"
multiqc "${filename}/rawdata/" -o "${filename}/reports/" -n "${filename}_multiqc_report_postqc"

echo "MultiQC report generated in ${filename}/reports/${filename}_multiqc_report_postqc.html"

echo "fastp processing completed for all SRA IDs."

echo "End of the script"
