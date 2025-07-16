#!/bin/bash

set -euo pipefail


inputdir="/media/labatl/HD-PCGU3-A/2506_processed"
outputdir="/home/labatl/devprojects/2507_breakend"

# Input parameters


nproc=16  # adjust based on available cores
min_mapq=30  # minimum mapping quality
htlv_ref="J02029.1"  # HTLV-1 reference name (adjust if using different HTLV)
java_mem="-Xmx24G -XX:MaxDirectMemorySize=32G"

# Reading BAM files folder ending with _sorted.bam

echo "✅ Reading FASTQ files from ${inputdir} and generating initial list..."
for file in ${inputdir}/*__sorted.bam; do
    echo "$file" 
done | sort -u > "${outputdir}/project_list.txt"
echo "✅ Initial list of unique BAM sample prefixes saved to ${outputdir}/project_list.txt"

while IFS= read -r sra; do
    newname=$(basename "$sra" | sed -E 's/_S[0-9]{1,3}_.*$//')
    echo "✅ New name for SRA ${sra} is ${newname}"

    
    if [ "${outputdir}/${newname}.bam" -ef "${sra}" ]; then
    echo "✅ Fastq files already copied for ${newname}. Skipping copy step"
    else
        echo "...Copying and renaming ${sra} to ${outputdir}/${newname}.bam"
        rsync -av --inplace --no-whole-file "$sra" "${outputdir}/${newname}.bam"
        rsync -av --inplace --no-whole-file "$sra.bai" "${outputdir}/${newname}.bam.bai"
        echo "✅ BAM file for ${newname} copied to ${outputdir}/"
    fi
    bam_input="${outputdir}/${newname}.bam"
    # 1. Extract pure HTLV-mapped reads (high MAPQ, sorted & indexed)
    echo "Extracting pure HTLV-mapped reads..."
    samtools view -@ ${nproc} -b -q ${min_mapq} "${bam_input}" "${htlv_ref}" | \
    samtools sort -@ ${nproc} -o "${outputdir}/${newname}_htlv_pure.bam"
    samtools index -@ ${nproc} "${outputdir}/${newname}_htlv_pure.bam"

    # 2. Extract chimeric human-HTLV reads
    echo "Extracting chimeric human-HTLV reads..."
    # A. Unmapped reads with mate mapped to HTLV
    samtools view -@ ${nproc} -f 4 -F 8 "${bam_input}" | \
    awk -v ref="${htlv_ref}" '$6 ~ /[HS]/ && $7 == ref {print $1}' | \
    sort -u > "${outputdir}/${newname}_chimeric_readnames.txt"

    # B. Reads mapped to HTLV with mate mapped to human
    samtools view -@ ${nproc} "${bam_input}" "${htlv_ref}" | \
    awk '$7 ~ /^(chr)?[0-9XYMT]+$/ {print $1}' | \
    sort -u >> "${outputdir}/${newname}_chimeric_readnames.txt"

    # C. Extract all chimeric reads (sorted & indexed)
    samtools view -@ ${nproc} -b -N "${outputdir}/${newname}_chimeric_readnames.txt" "${bam_input}" | \
    samtools sort -@ ${nproc} -o "${outputdir}/${newname}_htlv_chimeric.bam"
    samtools index -@ ${nproc} "${outputdir}/${newname}_htlv_chimeric.bam"

    # 3. Merge pure and chimeric BAMs (sorted & indexed)
    echo "Merging pure and chimeric reads..."
    samtools merge -@ ${nproc} -f "${outputdir}/${newname}_htlv_all.bam" \
    "${outputdir}/${newname}_htlv_pure.bam" \
    "${outputdir}/${newname}_htlv_chimeric.bam"
    samtools index -@ ${nproc} "${outputdir}/${newname}_htlv_all.bam"


    echo "✅ Processing for ${newname} completed."

    # 5. Mark duplicates
    echo "Marking duplicates..."
    gatk --java-options "${java_mem}" MarkDuplicates \
    -I "${outputdir}/${newname}_htlv_all.bam" \
    -O "${outputdir}/${newname}_htlv_dedup.bam"  \
    -M ${outputdir}/${newname}_dedup_metrics.txt \
    --CREATE_INDEX false \
    --REMOVE_DUPLICATES true

    samtools index -@ ${nproc} "${outputdir}/${newname}_htlv_dedup.bam"
    # 4. Cleanup
    rm -rf  "${outputdir}/${newname}_chimeric_readnames.txt" \
           "${outputdir}/${newname}_htlv_pure.bam" \
           "${outputdir}/${newname}_htlv_pure.bam.bai" \
            "${outputdir}/${newname}_htlv_chimeric.bam" \
            "${outputdir}/${newname}_htlv_chimeric.bam.bai" \
            "${bam_input}" \
            "${bam_input}.bai" \
            "${outputdir}/${newname}_htlv_all.bam" \
            "${outputdir}/${newname}_htlv_all.bam.bai" 

done < "${outputdir}/project_list.txt"







