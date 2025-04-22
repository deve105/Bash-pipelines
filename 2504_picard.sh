#!/bin/bash

set -euo pipefail #Stric mode


maindir="/home/labatl/devprojects"
platform="illumina"
indexgenome="/home/labatl/devapps/2407_references/2405project_hg38/2405_hg38htlv.fa"
proj_dir="/home/labatl/devprojects/Peru_IRID"
htlvgenome="J02029.1"
thrx=42

# It has to be executed in the folder where all the bam are, to create a list of files

for file in *.bam; do
    echo "$file" | sed 's/.bam$//'
done | uniq > "${maindir}/Peru_IRID/HTLV_1.txt"

echo "Initial list of renamed files in ${maindir}/Peru_IRID/HTLV_1.txt"


#### Looping the text file
while IFS= read -r sample_ID; do
    echo "Extracting HTLV-1 specific reads for ${sample_ID}" 
    # Specific Reads
    echo "Extracting HTLV-1-specific reads"
    samtools view -b -F 4 -@ ${thrx} \
        ${proj_dir}/bam/${sample_ID}.bam \
        ${htlvgenome} > \
        ${proj_dir}/htlv1/${sample_ID}_htlv1_specific.bam

    # Chimeric reads
    echo "Extracting discordant read pairs"
    samtools view -b -F 4 -f 0x1 -F 0x2 -@ ${thrx} \
        ${proj_dir}/bam/${sample_ID}.bam \
        ${htlvgenome} > \
        ${proj_dir}/htlv1/${sample_ID}_discordant.bam

    # Split reads
    echo "Extracting split reads"
    samtools view -h ${proj_dir}/bam/${sample_ID}.bam ${htlvgenome} | \
    awk '$6 ~ /[0-9]+[NMS]/ || $7 != "="' | \
    samtools view -b -@ ${thrx} > \
        ${proj_dir}/htlv1/${sample_ID}_split.bam

    # Merging Specific htlv reads, chimeric reads and split reads
    echo "Merging all HTLV-1 reads"
    samtools merge -f -@ ${thrx} \
        ${proj_dir}/htlv1/${sample_ID}_htlv1_all.bam \
        ${proj_dir}/htlv1/${sample_ID}_htlv1_specific.bam \
        ${proj_dir}/htlv1/${sample_ID}_discordant.bam \
        ${proj_dir}/htlv1/${sample_ID}_split.bam
    rm -rf ${proj_dir}/htlv1/${sample_ID}_htlv1_specific.bam \
        ${proj_dir}/htlv1/${sample_ID}_discordant.bam \
        ${proj_dir}/htlv1/${sample_ID}_split.bam
    # Sort and Index
    samtools sort -@ ${thrx} -m 2G\
        ${proj_dir}/htlv1/${sample_ID}_htlv1_all.bam \
        -o ${proj_dir}/htlv1/${sample_ID}_htlv1_final.sorted.bam
    samtools index ${proj_dir}/htlv1/${sample_ID}_htlv1_final.sorted.bam
    rm -rf ${proj_dir}/htlv1/${sample_ID}_htlv1_all.bam
    # Picard to mark duplicates
    picard MarkDuplicates \
    -I ${proj_dir}/htlv1/${sample_ID}_htlv1_final.sorted.bam \
    -O ${proj_dir}/htlv1/${sample_ID}_htlv1_marked.bam \
    -M ${proj_dir}/htlv1/${sample_ID}_htlv1_dup_metrics.txt \
    -CREATE_INDEX true \
    -VALIDATION_STRINGENCY SILENT \
    -REMOVE_DUPLICATES true
    rm -rf ${proj_dir}/htlv1/${sample_ID}_htlv1_final.sorted.bam
    # Coverage analysis with mosdepth
    echo "Running mosdepth for coverage analysis"
    mosdepth -t ${thrx} -m --by ${proj_dir}/J020209_1.fasta.bed \
        ${proj_dir}/mosdepth/${sample_ID}_v5_ \
        ${proj_dir}/htlv1/${sample_ID}_htlv1_marked.bam
        
    #${proj_dir}/J020209_1.fasta.bed -m --fast-mode 
    # Additional QC: Insert size metrics (for paired-end data)
    #picard CollectInsertSizeMetrics \
    #    I=${proj_dir}/htlv1/${sample_ID}_htlv1_marked.bam \
    #    O=${proj_dir}/htlv1/${sample_ID}_insert_size_metrics.txt \
    #    H=${proj_dir}/htlv1/${sample_ID}_insert_size_histogram.pdf \
    #    M=0.5
   
done < "${maindir}/Peru_IRID/HTLV_1.txt"

echo "End of the script"