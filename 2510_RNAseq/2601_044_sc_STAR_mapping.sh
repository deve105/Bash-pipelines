#!/bin/bash

set -euo pipefail #Stric mode

echo "Script for STAR mapping" 
echo "by Daniel Enriquez Version V1.2 2026-01-07"


hdd_tmp="/media/htlvatl/mem"
filename=$(basename "$1")
filename="${filename%.*}" # Removes the last extension
filename_path="${hdd_tmp}/${filename}"
numthreads=16
overhang=14

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

for tool in STAR samtools multiqc; do
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
        STAR --twopassMode Basic \
                --runThreadN "${numthreads}" \
                --genomeDir "${index_dir}" \
                --readFilesIn "${fastq1}" "${fastq2}" \
                --readFilesCommand zcat \
                --sjdbOverhang "${overhang}" \
                --outSAMtype BAM SortedByCoordinate \
                --sjdbGTFfile "${gtf_file}"  \
                --twopass1readsN -1 \
                --quantMode TranscriptomeSAM GeneCounts \
                --outSAMattrRGline ID:$sra SM:$sra PL:$platform \
                --outFileNamePrefix "${bamq}" \
                --outFilterMultimapScoreRange 1 \
                --outFilterMultimapNmax 20 \
                --outFilterMismatchNmax 10 \
                --alignIntronMax 500000 \
                --alignMatesGapMax 1000000 \
                --sjdbScore 2 \
                --alignSJDBoverhangMin 1 \
  	            --outFilterMatchNminOverLread 0.33 \
                --outFilterScoreMinOverLread 0.33 \
                --outSAMattributes NH HI NM MD AS XS nM \
                --outSAMunmapped Within \
                --chimSegmentMin 12 \
                --chimJunctionOverhangMin 8 \
                --chimOutJunctionFormat 1 \
                --alignSJstitchMismatchNmax 5 -1 5 5 \
                --chimMultimapScoreRange 3 \
                --chimScoreJunctionNonGTAG -4 \
                --chimMultimapNmax 20 \
                --chimNonchimScoreDropMin 10 \
                --peOverlapNbasesMin 12 \
                --peOverlapMMp 0.1 \
                --alignInsertionFlush Right \
                --alignSplicedMateMapLminOverLmate 0 \
                --alignSplicedMateMapLmin 30 \
             2>&1 | tee -a "${filename_path}/logs/${sra}_qc_star.log"
        samtools index ${SAMPLE}Aligned.sortedByCoord.out.bam
        rm -r ${SAMPLE}_STARtmp

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
