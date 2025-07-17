#!/bin/bash

set -euo pipefail



outputdir="/home/labatl/devprojects/2507_breakend"

# Input parameters


nproc=16  # adjust based on available cores
min_mapq=30  # minimum mapping quality
htlv_ref="/home/labatl/devapps/2407_references/2405project_hg38/J20209_1.fasta"  
# HTLV-1 reference name (adjust if using different HTLV)
ref_genome="/home/labatl/devapps/2407_references/2405project_hg38/2405_hg38_bwa/2405_hg38htlv.fa"
gridss_jar="/home/labatl/devapps/gridss/gridss.jar"
virusbreakenddb="/home/labatl/devapps/2407_references/2405project_hg38/gridss/virusbreakenddb_20210401"
# Reading BAM files folder ending with _sorted.bam


while IFS= read -r sra; do
    newname=$(basename "$sra" | sed -E 's/_S[0-9]{1,3}_.*$//')
    echo "✅ New name for SRA ${sra} is ${newname}"

    # virus bam aligned
    virusbam="${outputdir}/${newname}_htlv_dedup.bam"
    finaloutput="${outputdir}/output/${newname}_breakend.vcf"
    
    # virusbreakend with default virus database
    virusbreakend \
        -r "${ref_genome}" \
        -j "${gridss_jar}" \
        -o "${finaloutput}" \
        -w "${outputdir}/tmp" \
        --db "${virusbreakenddb}" \
        --threads 4 \
        --kraken2args "--memory-mapping" \
        "${virusbam}"
#--viralreferences "/home/labatl/devprojects/2507_breakend/htlv.txt" \
done < "${outputdir}/project_list.txt"