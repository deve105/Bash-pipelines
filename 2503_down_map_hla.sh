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

# Step 1 - create a folder with the project ID
filename=$(basename "$1")
filename="${filename%.*}" # Removes the last extension

echo "Step 1. Create a folder for the project '$filename'"
if [ -d "$filename" ]; then
	echo "Folder '$filename' already exists. Be careful"
	mkdir -p "${filename}"/{rawdata,bam,reports,hlala}
else
	echo "New folder '$filename' created."
	mkdir -p "${filename}"/{rawdata,bam,reports,hlala} #Creates the folder structure
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
    if ! fasterq-dump --mem 100G \
                      --threads 32 \
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
        gzip -f "${filename}/rawdata/${sra}_1.fq" "${filename}/rawdata/${sra}_2.fq"
    elif [ -f "${filename}/rawdata/${sra}.fastq" ]; then
        # Single-end reads
        echo "Single-end reads detected for $sra"
        mv "${filename}/rawdata/${sra}.fastq" "${filename}/rawdata/${sra}_1.fq"
        gzip -f "${filename}/rawdata/${sra}_1.fq"
    else
        echo "Error: No valid FASTQ files found for $sra"
        exit 1
    fi
    echo "${sra} was downloaded and processed :)"
done < "$1"

# Step 4 - Run fastp on the rawdata folder
echo "Running fastp on the rawdata folder"
while IFS= read -r sra; do
    if [ -f "${filename}/rawdata/${sra}_1.fq.gz" ] && [ -f "${filename}/rawdata/${sra}_2.fq.gz" ]; then
    # Paired-end mode
        fastp -i "${filename}/rawdata/${sra}_1.fq.gz" \
		    -I "${filename}/rawdata/${sra}_2.fq.gz" \
            -3 \
            -o "${filename}/rawdata/${sra}_postqc_1.fq.gz" \
            -O "${filename}/rawdata/${sra}_postqc_2.fq.gz" \
			-j "${filename}/reports/${sra}_fastp.json" \
            -h "${filename}/reports/${sra}_fastp.html" \
            -q 30 \
            -u 40 \
            -n 5 \
            --detect_adapter_for_pe \
            -w 42
        # Removing extra files
		rm "${filename}/rawdata/${sra}_1.fq.gz" "${filename}/rawdata/${sra}_2.fq.gz"
        # New References
        fastq1="${filename}/rawdata/${sra}_postqc_1.fq.gz"
        fastq2="${filename}/rawdata/${sra}_postqc_2.fq.gz"
        bamq="${filename}/bam/${sra}"
        indexgenome="/home/labatl/devapps/2407_references/1000g38/genome.fa"
        #MAPPING
        echo "BWA-MEM2 mapping for ${sra}}"
        bwa-mem2 mem -t 42 "${indexgenome}" "${fastq1}" "${fastq2}" | \
        samtools view -o "${bamq}.bam"
        #SORTING
        echo "Sorting BAM file"
        samtools sort -@ 42 -o "${bamq}_sorted.bam" "${bamq}.bam"
        #INDEXING
        echo "Indexing BAM file"
        samtools index "${bamq}_sorted.bam"
        #HLA Typing
        echo "HLA Typing"
        HLA-LA.pl --BAM "${bamq}_sorted.bam" \
        --graph PRG_MHC_GRCh38_withIMGT \
        --sampleID ${sra} \
        --maxThreads 42 \
        --workingDir "${filename}/hlala"
        # Removing extra data    
        rm -rf "${bamq}.bam" "${bamq}_sorted.bam" "${bamq}_sorted.bam.bai"
        rm -rf "${fastq1}" "${fastq2}" 
        rm -rf ${filename}/hlala/${sra}/*.bam 
        rm -rf ${filename}/hlala/${sra}/*.bam.bai 
        rm -rf ${filename}/hlala/${sra}/*.fastq
        rm -rf ${filename}/hlala/${sra}/*.txt

    elif [ -f "${filename}/rawdata/${sra}_1.fq.gz" ]; then
        # Single-end mode
        fastp -i "${filename}/rawdata/${sra}_1.fq.gz"  \
		    -o "${filename}/rawdata/${sra}_postqc_1.fq.gz" \
            -j "${filename}/reports/${sra}_fastp.json" \
			-h "${filename}/reports/${sra}_fastp.html" \
            -q 30 \
            -u 40 \
            -n 5 \
            --detect_adapter_for_pe \
            -w 42
		rm "${filename}/rawdata/${sra}_1.fq.gz"
        # New References
        fastq1="${filename}/rawdata/${sra}_postqc_1.fq.gz"
        bamq="${filename}/bam/${sra}"
        indexgenome="/home/labatl/devapps/2407_references/1000g38/genome.fa"
        #MAPPING
        echo "BWA-MEM2 mapping for ${sra}}"
        bwa-mem2 mem -t 42 "${indexgenome}" "${fastq1}" | \
        samtools view -o "${bamq}.bam"
        #SORTING
        echo "Sorting BAM file"
        samtools sort -@ 42 -o "${bamq}_sorted.bam" "${bamq}.bam"
        #INDEXING
        echo "Indexing BAM file"
        samtools index "${bamq}_sorted.bam"
        #HLA Typing
        echo "HLA Typing"
        HLA-LA.pl --BAM "${bamq}_sorted.bam" \
        --graph PRG_MHC_GRCh38_withIMGT \
        --sampleID ${sra} \
        --maxThreads 42 \
        --workingDir "${filename}/hlala"
        # Removing extra data    
        rm -rf "${bamq}.bam" "${bamq}_sorted.bam" "${bamq}_sorted.bam.bai"
        rm -rf "${fastq1}" "${fastq2}" 
        rm -rf ${filename}/hlala/${sra}/*.bam 
        rm -rf ${filename}/hlala/${sra}/*.bam.bai 
        rm -rf ${filename}/hlala/${sra}/*.fastq
        rm -rf ${filename}/hlala/${sra}/*.txt 
    else
        echo "Error: No valid FASTQ files found for $sra"
        exit 1
    fi
done < "$1"

# Step 6 - Re-Run MultiQC on the rawdata folder
echo "Running MultiQC on the rawdata folder"
multiqc "${filename}/reports/" -o "${filename}/reports/" -n "${filename}_multiqc_report_postqc"

echo "MultiQC report generated in ${filename}/reports/${filename}_multiqc_report_postqc.html"

echo "fastp processing completed for all SRA IDs."

echo "End of the script"
