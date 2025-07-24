#!/bin/bash

set -euo pipefail


inputdir="/media/labatl/HD-PCGU3-A/2506_processed"
outputdir="/home/labatl/devprojects/2507_breakendv5"

# Input parameters
mkdir -p "${outputdir}/tmp" 

nproc=16  # adjust based on available cores
min_mapq=30  # minimum mapping quality

java_mem="-Xmx24G -XX:MaxDirectMemorySize=32G"

htlv_ref="/home/labatl/devapps/2407_references/2405project_hg38/J20209_1.fasta"  
# HTLV-1 reference name (adjust if using different HTLV)
ref_genome="/home/labatl/devapps/2407_references/2405project_hg38/2405_hg38_bwa/2405_hg38htlv.fa"
gridss_jar="/home/labatl/devapps/gridss/gridss.jar"
virusbreakenddb="/home/labatl/Downloads/virusbreakenddb_20210401_new"
listvirus="/home/labatl/devprojects/2507_breakend/htlv.txt"

onedrive="/home/labatl/OneDrive"
# Reading BAM files folder ending with _sorted.bam

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
    
    # 5. Mark duplicates
    echo "Marking duplicates..."
    gatk --java-options "${java_mem}" MarkDuplicates \
    -I "${outputdir}/${newname}.bam"\
    -O "${outputdir}/${newname}_dedup.bam"  \
    -M ${outputdir}/${newname}_dedup_metrics.txt \
    --CREATE_INDEX false \
    --REMOVE_DUPLICATES false

    samtools index -@ ${nproc} "${outputdir}/${newname}_dedup.bam"

    # virus bam aligned
    virusbam="${outputdir}/${newname}_dedup.bam"
    finaloutput="${onedrive}/output2507_meta/${newname}_breakend.vcf"
    
    # Run virusbreakend, continue to next file if it fails
    if ! virusbreakend \
        -r "${ref_genome}" \
        -j "${gridss_jar}" \
        -o "${finaloutput}" \
        -w "${outputdir}/tmp" \
        --db "${virusbreakenddb}" \
        --threads 24 \
        --kraken2args "--memory-mapping" \
        --minviralcoverage 1 \
        --minreads 1 \
        --host human \
        --viralreferences "${listvirus}" \
        "${sra}"; then
        echo "❌ VirusBreakend failed for ${newname}, skipping to next."
    fi
    # 4. Cleanup
    rm -rf  "${bam_input}" \
            "${bam_input}.bai" \
            "${outputdir}/${newname}_dedup.bam" \
            "${outputdir}/${newname}_dedup.bam.bai" 

done < "${outputdir}/project_list.txt"

onedrive -s






