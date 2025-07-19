#!/bin/bash

set -euo pipefail



outputdir="/home/labatl/devprojects/2507_breakendv2"

# Input parameters


nproc=16  # adjust based on available cores
min_mapq=30  # minimum mapping quality
htlv_ref="/home/labatl/devapps/2407_references/2405project_hg38/J20209_1.fasta"  
# HTLV-1 reference name (adjust if using different HTLV)
ref_genome="/home/labatl/devapps/2407_references/2405project_hg38/2405_hg38_bwa/2405_hg38htlv.fa"
gridss_jar="/home/labatl/devapps/gridss/gridss.jar"
virusbreakenddb="/home/labatl/Downloads/virusbreakenddb_20210401_new"
listvirus="/home/labatl/devprojects/2507_breakend/htlv.txt"
# Reading BAM files folder ending with _sorted.bam


while IFS= read -r sra; do
    newname=$(basename "$sra" | sed -E 's/_S[0-9]{1,3}_.*$//')
    echo "✅ New name for SRA ${sra} is ${newname}"

    # virus bam aligned
    virusbam="${outputdir}/${newname}_htlv_dedup.bam"
    finaloutput="${outputdir}/output4/${newname}_breakend.vcf"
    
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
        "${virusbam}"; then
        echo "❌ VirusBreakend failed for ${newname}, skipping to next."
    fi
#--viralreferences "/home/labatl/devprojects/2507_breakend/htlv.txt" \
done < "${outputdir}/project_list.txt"