#!/bin/bash

#set -euo pipefail

echo "========================================"
echo "✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅"
echo "Script for FASTQ virus and mutations by Daniel Enriquez-Vera"
echo "Version V1 2025-06-3"

## Input
inputdir="/home/htlvatl/Documents/devprojects/Peru_IRID/fastq"
outputdir="/home/htlvatl/Documents/devprojects/Peru_IRID"
ref_genome="/home/htlvatl/Documents/devapps/references/GRCh38_HTLV1.fa"
conda_env="devpipe"

## Folder structure inside outputdir

mkdir -p "${outputdir}/Peru_IRID"/{temp,fastp,bam,reports}

## Software paths *mamba environment fastp
bwamem2="/home/htlvatl/Documents/devapps/bwa-mem2-2.3_x64-linux/bwa-mem2.avx2"

########################### Software checks ###########################

### Check if BWA-MEM2 is installed and executable
if [ ! -f "$bwamem2" ]; then
    printf "❌ Error: BWA-MEM2 executable not found at: %s\n" "$bwamem2"
    echo "Please check if the path is correct or install BWA-MEM2"
    exit 1
elif [ ! -x "$bwamem2" ]; then
    printf "❌ Error: BWA-MEM2 is not executable: %s\n" "$bwamem2"
    echo "Run: chmod +x $bwamem2"
    exit 1
else
    printf "✅ BWA-MEM2: Found at %s\n" "$bwamem2"
    # Test if it actually works (improved testing)
    if "$bwamem2" 2>&1 | grep -qi "usage\|bwa-mem2\|version" || \
       "$bwamem2" mem 2>&1 | grep -qi "usage\|required\|options"; then
        printf "✅ BWA-MEM2: Executable and working\n"
    else
        printf "⚠️  Warning: BWA-MEM2 executable found but may not work properly\n"
    fi
fi

### Check if conda/mamba environment is activated
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda/mamba is not installed or not in PATH."
    echo "Please install Miniforge or Anaconda and activate your environment."
    exit 1
else
    echo "✅ Conda/mamba is available."
    if ! conda info --envs | grep -q "$conda_env"; then
        echo "❌ Error: '$conda_env' environment not found."
        echo "Please create it with: mamba create -n $conda_env -c bioconda -c conda-forge bwa-mem2 samtools fastp multiqc python=3.9 -y"
        exit 1
    else
        echo "✅ '$conda_env' environment is available."
    fi
fi  

### Check if the environment is currently activated
current_env=$(conda info --json | grep -o '"active_prefix_name": "[^"]*"' | cut -d'"' -f4)
if [ "$current_env" = "$conda_env" ]; then
    echo "✅ '$conda_env' environment is currently ACTIVATED."
else
    echo "⚠️  WARNING: '$conda_env' environment exists but is NOT activated."
    exit 1
fi  

### Check if fastp is installed
if ! command -v fastp &> /dev/null; then
    echo "❌ Error: fastp is not installed or not in PATH."
    echo "Please install it in your conda/mamba environment."
    exit 1
else
    echo "✅ fastp is available."
fi

### Check if multiqc is installed
if ! command -v multiqc &> /dev/null; then
    echo "❌ Error: multiqc is not installed or not in PATH."
    echo "Please install it in your conda/mamba environment."
    exit 1
else
    echo "✅ multiqc is available."
fi

### Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "❌ Error: samtools is not installed or not in PATH."
    echo "Please install it in your conda/mamba environment."
    exit 1
else
    echo "✅ samtools is available."
fi
#######################################################

# Generate a list text of all fastq files in the input directory
echo "✅ Reading FASTQ files from ${inputdir} and generating initial list..."
for file in ${inputdir}/*.fastq.gz; do
    echo "$file" | sed 's/R[0-9]_001.fastq.gz$//'
done | uniq > "${outputdir}/Peru_IRID_list.txt"
echo "✅ Initial list of renamed files in ${outputdir}/Peru_IRID_list.txt"

# Loop through the list of SRA files and process them

while IFS= read -r sra; do
    if [ ! -f "${inputdir}/${sra}R1_001.fastq.gz" ] || [ ! -f "${inputdir}/${sra}R3_001.fastq.gz" ]; then
        echo "❌ Error: Missing FASTQ files for ${sra}"
        echo "Expected ${sra}R1_001.fastq.gz and ${sra}R3_001.fastq.gz in ${inputdir}."
        exit 1
    fi
    echo "✅ Processing SRA: ${sra}"
    
    # Renaming the SRA file to a new name
    newname=$(echo "${sra}" | sed -E 's/^.*Nakahata_//; s/_S[0-9].*$//' | sed -E 's/^.*IRID/IRID/')
    echo "✅ New name for SRA ${sra} is ${newname}"
    
    # Copy and rename the fastq files
    echo "Copying and renaming ${sra}R1_001.fastq.gz to ${newname}_1.fq.gz"
    rsync -av "${inputdir}/${sra}R1_001.fastq.gz" "${outputdir}/Peru_IRID/temp/${newname}_1.fq.gz"
    echo "Copying and renaming ${sra}R3_001.fastq.gz to ${newname}_2.fq.gz"
    rsync -av "${inputdir}/${sra}R3_001.fastq.gz" "${outputdir}/Peru_IRID/temp/${newname}_2.fq.gz"
    echo "${newname}" >> "${outputdir}/Peru_IRID/Peru_copied.txt"
    echo "✅ Fastq_1 and Fastq_2 of ${newname} copied to ${outputdir}/Peru_IRID/temp/"

    # Quality control with fastp
    echo "Quality control for ${newname}"
    fastp -i "${outputdir}/Peru_IRID/temp/${newname}_1.fq.gz" \
        -I "${outputdir}/Peru_IRID/temp/${newname}_2.fq.gz" \
        -o "${outputdir}/Peru_IRID/temp/${newname}_postqc_1.fq.gz" \
        -O "${outputdir}/Peru_IRID/temp/${newname}_postqc_2.fq.gz" \
        -j "${outputdir}/Peru_IRID/fastp/${newname}_fastp.json" \
        -h "${outputdir}/Peru_IRID/fastp/${newname}_fastp.html" \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 40 \
        --n_base_limit 5 \
        --length_required 50 \
        --cut_right \
        --correction \
        --cut_mean_quality 20 \
        --cut_window_size 4 \
        --detect_adapter_for_pe \
        --dont_eval_duplication \
        --thread ${nproc}
    echo "✅ Quality control completed for ${newname}"
	
    rm "${outputdir}/Peru_IRID/temp/${newname}_1.fq.gz" \
        "${outputdir}/Peru_IRID/temp/${newname}_2.fq.gz"
    echo "✅ ${newname}_1.fq.gz and ${newname}_2.fq.gz removed from temp folder"

    fastq1="${outputdir}/Peru_IRID/temp/${newname}_postqc_1.fq.gz"
    fastq2="${outputdir}/Peru_IRID/temp/${newname}_postqc_2.fq.gz"
    bamq="${outputdir}/Peru_IRID/bam/${newname}"
    
    #MAPPING
    echo "BWA-MEM2 mapping for ${newname}"
    "$bwamem2" mem \
        -t ${nproc} \
        -R "@RG\tID:${newname}\tSM:${newname}\tPL:ILLUMINA" \
        "${ref_genome}" \
        "${fastq1}" \
        "${fastq2}" | \
    samtools view -@ 8 -bS - > "${bamq}.bam"
    echo "✅ BAM file created: ${bamq}.bam"

    #SORTING
    echo "Sorting BAM file"
    samtools sort -@ ${nproc} -o "${bamq}_sorted.bam" "${bamq}.bam"
    echo "✅ Sorted BAM file created: ${bamq}_sorted.bam"

    #INDEXING
    echo "Indexing BAM file"
    samtools index "${bamq}_sorted.bam"
    echo "✅ Indexed BAM file created: ${bamq}_sorted.bam.bai"

    rm -f "${fastq1}" "${fastq2}" "${bamq}.bam"
    echo "✅ Intermediate BAM file cleaned up"

    # Picard to mark duplicates
    picard MarkDuplicates \
    -I ${proj_dir}/htlv1/${sample_ID}_htlv1_final.sorted.bam \
    -O ${proj_dir}/htlv1/${sample_ID}_htlv1_marked.bam \
    -M ${proj_dir}/htlv1/${sample_ID}_htlv1_dup_metrics.txt \
    -CREATE_INDEX true \
    -VALIDATION_STRINGENCY SILENT \
    -REMOVE_DUPLICATES true
    
    rm -rf ${proj_dir}/htlv1/${sample_ID}_htlv1_final.sorted.bam


done < "${maindir}/Peru_IRID/Peru_IRID_list.txt"



# FIX MULTIQC:
echo "📊 Running MultiQC"
multiqc "${outputdir}/Peru_IRID/fastp" -o "${outputdir}/Peru_IRID/reports" -n "multiqc_report_postqc"
echo "✅ MultiQC report: ${outputdir}/Peru_IRID/multiqc_report_postqc.html"