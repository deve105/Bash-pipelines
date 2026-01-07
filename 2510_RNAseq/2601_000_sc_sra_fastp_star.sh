#!/bin/bash

set -euo pipefail #Stric mode

echo "Script for SRA download, fastp QC, STAR mapping, and multiQC" 
echo "by Daniel Enriquez Version V1.2 2026-01-07"


hdd_tmp="/media/htlvatl/mem"
filename=$(basename "$1")
filename="${filename%.*}" # Removes the last extension
filename_path="${hdd_tmp}/${filename}"

# STAR index directory
ref_dir="/home/htlvatl/coreapps/2025_refgenomes/2601_001_di_htlv1_ebv_hg38"
merged_fasta="${ref_dir}/2601_001_rf_htlv1_ebv_hg38.fa"
merged_gtf="${ref_dir}/2601_001_rf_htlv1_ebv_hg38.gtf"
index_dir="${ref_dir}/2601_001_rf_star_index_merged"

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

# Step 1 create project folder
echo "Step 1. Create a folder for the project '$filename'"
if [ -d "$filename_path" ]; then
	echo "Folder '$filename_path' already exists. Be careful"
	mkdir -p "${filename_path}"/{rawdata,bam,logs,reports}
else
	echo "New folder '$filename_path' created."
	mkdir -p "${filename_path}"/{rawdata,bam,logs,reports} #Creates the folder structure
fi


if [ ! -d "$index_dir" ]; then
    printf "❌ Error: STAR index not found: %s\n" "$index_dir"
    printf "   Please provide STAR index directory as second argument\n"
    exit 1
fi

# ============================================================================
# STEP 1 - CHECK DEPENDENCIES
# ============================================================================

echo "🔍 Checking dependencies..."

for tool in prefetch fasterq-dump fastp STAR samtools multiqc; do
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
    cd "${hdd_tmp}"
    echo "Downloading ${sra}"

    # Step 3.1 - Prefetch the SRA data
    if ! prefetch "${sra}"; then
        echo "Error: prefetch failed for ${sra}"
        exit 1
    fi

    # Step 3.2 - Convert SRA to FASTQ using fasterq-dump
    if ! fasterq-dump --mem 30G \
                      --threads 12 \
                      --skip-technical \
                      --outdir "${filename_path}/rawdata/" \
                      --temp "${filename_path}" \
                      --split-3 \
                      "${sra}"; then
        echo "Error: fasterq-dump failed for ${sra}"
        exit 1
    fi

		
	# Check if the files are paired-end or single-end
    if [ -f "${filename_path}/rawdata/${sra}_1.fastq" ] && [ -f "${filename_path}/rawdata/${sra}_2.fastq" ]; then
        # Paired-end reads
        echo "Paired-end reads detected for $sra"
        mv "${filename_path}/rawdata/${sra}_1.fastq" "${filename_path}/rawdata/${sra}_1.fq"
        mv "${filename_path}/rawdata/${sra}_2.fastq" "${filename_path}/rawdata/${sra}_2.fq"
        gzip -f "${filename_path}/rawdata/${sra}_1.fq" "${filename_path}/rawdata/${sra}_2.fq"
    elif [ -f "${filename_path}/rawdata/${sra}.fastq" ]; then
        # Single-end reads
        echo "Single-end reads detected for $sra"
        mv "${filename_path}/rawdata/${sra}.fastq" "${filename_path}/rawdata/${sra}_1.fq"
        gzip -f "${filename_path}/rawdata/${sra}_1.fq"
    else
        echo "Error: No valid FASTQ files found for $sra"
        exit 1
    fi
    echo "${sra} was downloaded and processed :)"
    if [ -f "${filename_path}/rawdata/${sra}_1.fq.gz" ] && [ -f "${filename_path}/rawdata/${sra}_2.fq.gz" ]; then
        echo "Performing fastp quality control for paired-end reads: ${sra}"
    # Paired-end mode
        fastp -i "${filename_path}/rawdata/${sra}_1.fq.gz" \
		    -I "${filename_path}/rawdata/${sra}_2.fq.gz" \
            -o "${filename_path}/rawdata/${sra}_postqc_1.fq.gz" \
            -O "${filename_path}/rawdata/${sra}_postqc_2.fq.gz" \
			-j "${filename_path}/reports/${sra}_fastp.json" \
            -h "${filename_path}/reports/${sra}_fastp.html" \
            --detect_adapter_for_pe \
            --length_required 15 \
            --correction \
            -e 15 \
            -q 15 \
            -w 12 \
            -u 40 \
            -n 5 \
            2>&1 | tee -a "${filename_path}/logs/${sra}_fastp.log"

    elif [ -f "${filename_path}/rawdata/${sra}_1.fq.gz" ]; then
        echo "Performing fastp quality control for single-end reads: ${sra}"
        # Single-end mode
        fastp -i "${filename_path}/rawdata/${sra}_1.fq.gz" \
            -o "${filename_path}/rawdata/${sra}_postqc_1.fq.gz" \
			-j "${filename_path}/reports/${sra}_fastp.json" \
            -h "${filename_path}/reports/${sra}_fastp.html" \
            --length_required 15 \
            --correction \
            -e 15 \
            -q 15 \
            -w 12 \
            -u 40 \
            -n 5 \
            2>&1 | tee -a "${filename_path}/logs/${sra}_fastp.log"
    else
        echo "Error: No valid FASTQ files found for $sra"
        exit 1
    fi
    if [ -f "${filename_path}/rawdata/${sra}_postqc_1.fq.gz" ] && [ -f "${filename_path}/rawdata/${sra}_postqc_2.fq.gz" ]; then
        echo "Mapping for paired-end reads: ${sra}"
        fastq1="${filename_path}/rawdata/${sra}_postqc_1.fq.gz"
        fastq2="${filename_path}/rawdata/${sra}_postqc_2.fq.gz"
        fastr1="${filename_path}/rawdata/${sra}_1.fq.gz"
        fastr2="${filename_path}/rawdata/${sra}_2.fq.gz"
        bamq="${filename_path}/bam/${sra}_postqc"
        bamraw="${filename_path}/bam/${sra}_raw"
        #MAPPING
        echo "STAR mapping for ${sra} for post-QC reads"
        # Paired-end alignment
        STAR --runThreadN 12 \
             --genomeDir "${index_dir}" \
             --readFilesIn "${fastq1}" "${fastq2}" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${bamq}" \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard \
             --sjdbOverhang 99 \
             --quantMode TranscriptomeSAM GeneCounts \
             --twopassMode Basic \
             2>&1 | tee -a "${filename_path}/logs/${sra}_qc_star.log"

        echo "STAR mapping for ${sra} for raw reads"
        STAR --runThreadN 12 \
             --genomeDir "${index_dir}" \
             --readFilesIn "${fastr1}" "${fastr2}" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${bamraw}" \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard \
             --sjdbOverhang 99 \
             --quantMode TranscriptomeSAM GeneCounts \
             --twopassMode Basic \
             2>&1 | tee -a "${filename_path}/logs/${sra}_raw_star.log"
    elif [ -f "${filename_path}/rawdata/${sra}_postqc_1.fq.gz" ]; then 
    echo "Mapping for single-end reads: ${sra}"
        fastq1="${filename_path}/rawdata/${sra}_postqc_1.fq.gz"
        fastr1="${filename_path}/rawdata/${sra}_1.fq.gz"
        bamq="${filename_path}/bam/${sra}_postqc"
        bamraw="${filename_path}/bam/${sra}_raw"
    echo "STAR mapping for ${sra} for post-QC reads"
        # Paired-end alignment
        STAR --runThreadN 12 \
             --genomeDir "${index_dir}" \
             --readFilesIn "${fastq1}" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${bamq}" \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard \
             --sjdbOverhang 99 \
             --quantMode TranscriptomeSAM GeneCounts \
             --twopassMode Basic \
             2>&1 | tee -a "${filename_path}/logs/${sra}_qc_star.log"

        echo "STAR mapping for ${sra} for raw reads"
        STAR --runThreadN 12 \
             --genomeDir "${index_dir}" \
             --readFilesIn "${fastr1}" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${bamraw}" \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard \
             --sjdbOverhang 99 \
             --quantMode TranscriptomeSAM GeneCounts \
             --twopassMode Basic \
             2>&1 | tee -a "${filename_path}/logs/${sra}_raw_star.log"
 
    fi

    ## Deleting intermediate files to save space
    rm -f "${filename_path}/rawdata/${sra}_1.fq.gz"
    rm -f "${filename_path}/rawdata/${sra}_2.fq.gz"
    rm -f "${filename_path}/rawdata/${sra}_postqc_1.fq.gz"
    rm -f "${filename_path}/rawdata/${sra}_postqc_2.fq.gz"

done < "$1"


