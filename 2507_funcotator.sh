#!/bin/bash

set -euo pipefail

echo "========================================"
echo "✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅"
echo "Script for FASTQ virus and mutations by Daniel Enriquez-Vera"
echo "Version V1 2025-06-3"

## Input
inputdir="/media/labatl/HD-PCGU3-A/2506_fastq"
outputdir="/home/labatl/devprojects/2506_fastq"
outdir2="/media/labatl/HD-PCGU3-A/2506_processed"
ref_genome="/home/labatl/devapps/2407_references/2405project_hg38/2405_hg38htlv.fa"
htlv1_ref="/home/labatl/devapps/2407_references/htlv1/J20209_1.fasta"
pon="/home/labatl/devapps/2407_references/1000g_pon.hg38.vcf.gz"
germline_resource="/home/labatl/devapps/2407_references/af-only-gnomad.hg38.vcf.gz"
funcotator_data="/home/labatl/devapps/2407_references/funcotator_dataSources.v1.8.hg38.20230908s"
known_sites="/home/labatl/devapps/2407_references/Homo_sapiens_assembly38.dbsnp138.vcf"

outputdir_req="/home/labatl/devprojects/2506_fastq/vcf_recalibration"
nproc=32
picard_cmd="/home/labatl/devapps/picard.jar"
gatk_cmd="/home/labatl/devapps/gatk/gatk"
java_mem="-Xmx24G -XX:MaxDirectMemorySize=32G"    

omni="/home/labatl/devapps/2407_references/1000G_omni2.5.hg38.vcf.gz"
dbsnp="/home/labatl/devapps/2407_references/funcotator_dataSources.v1.8.hg38.20230908s/dbsnp/hg38/hg38_All_20180418.vcf.gz"
thousandG="/home/labatl/devapps/2407_references/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
hapmap="/home/labatl/devapps/2407_references/hapmap_3.3.hg38.vcf.gz"

while IFS= read -r newname; do
    if false; then
    # This seems to only work after haplotype
    gatk VariantRecalibrator \
            -R "${ref_genome}" \
            -V "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz" \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousandG} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode BOTH \
            -O "${outputdir_req}/${newname}_output.tranches" \
            --tranches-file "${outputdir_req}/${newname}_output.tranches"

    gatk ApplyVQSR \
        -R "${ref_genome}" \
        -V "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz" \
        --recal-file "${outputdir_req}/${newname}_output.recal" \
        --tranches-file "${outputdir_req}/${newname}_output.tranches" \
        -mode BOTH \
        -ts-filter-level 95.0 \
        -O "${outputdir_req}/${newname}_recalibrated_variant.vcf"
    
    gatk VariantFiltration \
            -R ${ref_genome} \
            -V "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz"  \
            -O "${outputdir_req}/${newname}_final_filtered_variant.vcf"\
            -filter-name "QD_filter" -filter "QD < 5.0" \
            -filter-name "FS_filter" -filter "FS > 30.0" \
            -filter-name "MQ_filter" -filter "MQ < 30.0" \
            -filter-name "SOR_filter" -filter "SOR > 4.0" \
            -filter-name "MQRankSum_filter" -filter "MQRankSum < -5" \
            -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -5.0" \
            --filter-name "DP_filter" --filter "DP < 30.0" \
            --filter-name "AF_filter" --filter "AF < 0.05"    
    fi
    gatk VariantFiltration \
        -R "${ref_genome}" \
        -V "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz" \
        -O "${outputdir_req}/${newname}_final_filtered_variant.vcf" \
        --filter-name "low_AF" --filter-expression "vc.getAttribute('AF') < 0.01" \
        --filter-name "low_AD" --filter-expression "vc.getGenotype('TUMOR').getAD().1 < 5" \
        --filter-name "low_DP" --filter-expression "DP < 20" \
        --filter-name "strand_bias" --filter-expression "FS > 25.0" \
        --filter-name "weak_evidence" --filter-expression "TLOD < 3.0" \
        --filter-name "contamination" --filter-expression "CONTQ < 20" \
        --filter-name "germline_risk" --filter-expression "GERMQ < 20" \
        --filter-name "duplicate_artifact" --filter-expression "AS_UNIQ_ALT_READ_COUNT[0] < 2" \
        --filter-name "base_qual" --filter-expression "MBQ[0] < 25 || MBQ[1] < 25" \
        --filter-name "fragment_bias" --filter-expression "MFRL[1] - MFRL[0] > 100"  
    
    gatk --java-options "${java_mem}" Funcotator \
     --variant "${outputdir_req}/${newname}_final_filtered_variant.vcf" \
     --reference "${ref_genome}" \
     --output "${outputdir_req}/${newname}_final.maf" \
     --output-file-format MAF \
     --data-sources-path ${funcotator_data} \
     --ref-version hg38 \
     --transcript-selection-mode BEST_EFFECT \
     --remove-filtered-variants 

done < "/home/labatl/devprojects/2506_fastq/250711_foranalysis.txt"
