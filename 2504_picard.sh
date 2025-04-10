#!/bin/bash

set -euo pipefail #Stric mode


maindir="/home/labatl/devprojects"
platform="illumina"
indexgenome="/home/labatl/devapps/2407_references/2405project_hg38/2405_hg38htlv.fa"
proj_dir="/home/labatl/devprojects/Peru_IRID"
htlvgenome="J02029.1"
thrx=42

for file in *.bam; do
    echo "$file" | sed 's/.bam$//'
done | uniq > "${maindir}/Peru_IRID/HTLV_1.txt"

echo "Initial list of renamed files in ${maindir}/Peru_IRID/HTLV_1.txt.txt"

while IFS= read -r sra; do
    echo "Extracting HTLV-1 specific reads for ${sample_ID}" 
    samtools view -b -@ ${thrx} \
        ${proj_dir}/bam/${sample_ID}.bam \
        ${htlvgenome} > \
        ${proj_dir}/htlv1/${sample_ID}_htlv1.bam 
    picard MarkDuplicates \
        I=${proj_dir}/htlv1/${sample_ID}_htlv1.bam \
        O=${proj_dir}/htlv1/${sample_ID}_htlv1_marked.bam \
        M=${proj_dir}/htlv1/${sample_ID}_htlv1_dup_metrics.txt \
        CREATE_INDEX=true
    mosdepth -n --by htlv1.gff3 ${sample}.cov ${sample}.htlv1.marked.bam


    if [ ! -f "${sra}R1_001.fastq.gz" ] || [ ! -f "${sra}R3_001.fastq.gz" ]; then
        echo "Error: Missing FASTQ files for ${sra}. Expected ${sra}R1_001.fastq.gz and ${sra}R3_001.fastq.gz."
        exit 1
    fi
done < "${maindir}/Peru_IRID/HTLV_1.txt.txt"
echo "End of the script"