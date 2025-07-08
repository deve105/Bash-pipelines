#!/bin/bash

set -euo pipefail

echo "========================================"
echo "✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅"
echo "Script for FASTQ virus and mutations by Daniel Enriquez-Vera"
echo "Version V1 2025-06-3"

## Input
inputdir="/media/labatl/HD-PCGU3-A/test"
outputdir="/home/labatl/devprojects/test"
ref_genome="/home/labatl/devapps/2407_references/2405project_hg38/2405_hg38htlv.fa"
htlv1_ref="/home/labatl/devapps/2407_references/htlv1/J20209_1.fasta"
pon="/home/labatl/devapps/2407_references/1000g_pon.hg38.vcf.gz"
germline_resource="/home/labatl/devapps/2407_references/af-only-gnomad.hg38.vcf.gz"
funcotator_data="/home/labatl/devapps/2407_references/funcotator_dataSources.v1.7.20200521g"

nproc=32
picard_cmd="/home/labatl/devapps/picard.jar"
gatk_cmd="/home/labatl/devapps/gatk/gatk"
java_mem="-Xmx32G -XX:MaxDirectMemorySize=32G"

####Directories checks#########################
###############################################
### Check if input directory exists
if [ ! -d "$inputdir" ]; then
    echo "❌ Error: Input directory does not exist: $inputdir"
    exit 1
else
    echo "✅ Input directory exists: $inputdir"
fi

### Check if reference genome exists
#if [ ! -f "$ref_genome" ]; then
#    echo "❌ Error: Reference genome file does not exist: $ref_genome"
#    exit 1
#else
#    echo "✅ Reference genome file exists: $ref_genome"
#fi

### Check if output directory exists, create if not
if [ ! -d "$outputdir" ]; then
    echo "⚠️ Output directory does not exist. Creating: $outputdir"
    mkdir -p "$outputdir"
    if [ $? -ne 0 ]; then
        echo "❌ Error: Failed to create output directory: $outputdir"
        exit 1
    fi
else
    echo "✅ Output directory exists: $outputdir"
fi

## Folder structure inside outputdir
mkdir -p "${outputdir}/"{temp,fastp,bam,reports,virus,human,virus_variants,human_variants}

####Software checks############################
###############################################
### Check if BWA-MEM2 is installed and executable
for cmd in bwa-mem2 samtools fastp multiqc java; do
    if ! command -v $cmd &> /dev/null; then
        echo "❌ Error: $cmd is not installed or not in PATH"
        exit 1
    else
        echo "✅ $cmd is available"
    fi
done
###############################################
####Initial List###############################
###############################################
# Generate a list text of all fastq files in the input directory
date_suffix=$(date +"%y%m%d")

echo "✅ Reading FASTQ files from ${inputdir} and generating initial list..."
for file in ${inputdir}/*.fastq.gz; do
    echo "$file" | sed 's/R[0-9]_001.fastq.gz$//'
done | uniq > "${outputdir}/${date_suffix}_project_list.txt"

echo "✅ Initial list of renamed files in ${outputdir}/${date_suffix}_project_list.txt"

project_list="${outputdir}/${date_suffix}_project_list.txt"

###############################################
####Initial List###############################
###############################################
# Loop through the list of SRA files and process them

while IFS= read -r sra; do
    if [ ! -f "${sra}R1_001.fastq.gz" ] || [ ! -f "${sra}R3_001.fastq.gz" ]; then
        echo "❌ Error: Missing FASTQ files for ${sra}"
        echo "Expected ${sra}R1_001.fastq.gz and ${sra}R3_001.fastq.gz in ${inputdir}."
        exit 1
    fi
    echo "✅ Processing SRA: ${sra}"
    
    # Renaming the SRA file to a new name
    newname=$(echo "${sra}" | sed -E 's/^.*Nakahata_//; s/_S[0-9].*$//' | sed -E 's/^.*IRID/IRID/')
    echo "✅ New name for SRA ${sra} is ${newname}"
    
    # Copy and rename the fastq files
    echo "...Copying and renaming ${sra}R1_001.fastq.gz to ${newname}_1.fq.gz"
    rsync -av "${sra}R1_001.fastq.gz" "${outputdir}/temp/${newname}_1.fq.gz"
    echo "...Copying and renaming ${sra}R3_001.fastq.gz to ${newname}_2.fq.gz"
    rsync -av "${sra}R3_001.fastq.gz" "${outputdir}/temp/${newname}_2.fq.gz"
    echo "${newname}" >> "${outputdir}/${date_suffix}_foranalysis.txt"
    echo "${sra},${newname}" >> "${outputdir}/${date_suffix}_ref.txt"
    echo "✅ Fastq_1 and Fastq_2 of ${newname} copied to ${outputdir}/temp/"

    # Quality control with fastp
    echo "Quality control for ${newname}"
    fastp -i "${outputdir}/temp/${newname}_1.fq.gz" \
        -I "${outputdir}/temp/${newname}_2.fq.gz" \
        -o "${outputdir}/temp/${newname}_postqc_1.fq.gz" \
        -O "${outputdir}/temp/${newname}_postqc_2.fq.gz" \
        -j "${outputdir}/fastp/${newname}_fastp.json" \
        -h "${outputdir}/fastp/${newname}_fastp.html" \
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
	
    rm "${outputdir}/temp/${newname}_1.fq.gz" \
        "${outputdir}/temp/${newname}_2.fq.gz"
    echo "✅ ${newname}_1.fq.gz and ${newname}_2.fq.gz removed from temp folder"

    fastq1="${outputdir}/temp/${newname}_postqc_1.fq.gz"
    fastq2="${outputdir}/temp/${newname}_postqc_2.fq.gz"
    bamq="${outputdir}/bam/${newname}"
    
    #MAPPING
    echo "BWA-MEM2 mapping for ${newname}"
    bwa-mem2 mem \
        -t ${nproc} \
        -M \
        -R "@RG\tID:${newname}\tSM:${newname}\tPL:ILLUMINA" \
        "${ref_genome}" \
        "${fastq1}" \
        "${fastq2}" | \
    samtools view -@ ${nproc} -bS - > "${bamq}.bam"
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

    ###############################################
    ###############################################
    ###############################################
    
    echo "======== EXTRACTING HUMAN READS OF INTEREST ========"

    # Extract well-mapped human reads (properly paired, MAPQ>=30)
    echo "Extracting high-quality human reads..."
    # Filter reads mapped only to human chromosomes
    # Using regions directly in samtools (no temp file needed)

    # Create a temporary BED file with human chromosomes
    samtools idxstats "${bamq}_sorted.bam" | grep -v "J02029.1" | awk '{print $1"\t0\t"$2}' > "${outputdir}/human/${newname}_human_regions.bed"

    # Use the BED file to extract human reads
    samtools view -@ ${nproc} -b -f 2 -q 30 -M -L "${outputdir}/human/${newname}_human_regions.bed" \
    "${bamq}_sorted.bam" > "${outputdir}/human/${newname}_human.bam"

    # Sort and index
    samtools sort -@ ${nproc} -o "${outputdir}/human/${newname}_final.bam" "${outputdir}/human/${newname}_human.bam"
    samtools index "${outputdir}/human/${newname}_final.bam"
    echo "✅ Human reads extracted and sorted: ${outputdir}/human/${newname}_final.bam"
    
    # Mark duplicates
    echo "Marking duplicates in human reads..."
    java ${java_mem} -jar ${picard_cmd} MarkDuplicates \
        -I "${outputdir}/human/${newname}_final.bam" \
        -O "${outputdir}/human/${newname}_final_dedup.bam" \
        -M "${outputdir}/human/${newname}_dedup_metrics.txt" \
        -CREATE_INDEX true \
        -REMOVE_DUPLICATES true \
        -VALIDATION_STRINGENCY LENIENT

    
    echo "✅ Human reads processed and duplicates marked: ${outputdir}/human/${newname}_final_dedup.bam"
    rm -f "${outputdir}/human/${newname}_final.bam" \
          "${outputdir}/human/${newname}_final.bam.bai" \
          "${outputdir}/human/${newname}_human.bam" \
          "${outputdir}/human/${newname}_human_regions.bed"
    
    echo "✅ Cleaned up intermediate human BAM files"

    # Generate coverage and statistics
    echo "Generating statistics..."
    samtools flagstat "${outputdir}/human/${newname}_final_dedup.bam"  > \
        "${outputdir}/reports/${newname}.human.stats.txt"
    samtools coverage "${outputdir}/human/${newname}_final_dedup.bam"  > \
        "${outputdir}/reports/${newname}.human.coverage.txt"
    echo "✅ Statistics generated for human reads: ${outputdir}/reports/${newname}.human.stats.txt"

    #######################################

    # Run Mutect2 with panel of normals
    echo "Running Mutect2 for ${newname} with panel of normals..."
    gatk --java-options "${java_mem}" Mutect2 \
        -R "${ref_genome}" \
        -I "${outputdir}/human/${newname}_final_dedup.bam" \
        --panel-of-normals "${pon}" \
        --germline-resource "${germline_resource}" \
        --f1r2-tar-gz "${outputdir}/human_variants/${newname}_f1r2.tar.gz" \
        -O "${outputdir}/human_variants/${newname}_somatic_unfiltered.vcf.gz" \
        --max-population-af 0.01 \
        --genotype-germline-sites true \
        --genotype-pon-sites true
    
    # Learn read orientation model
    echo "Learning read orientation model for ${newname}..."
    gatk --java-options "${java_mem}" LearnReadOrientationModel \
        -I "${outputdir}/human_variants/${newname}_f1r2.tar.gz" \
        -O "${outputdir}/human_variants/${newname}_read_orientation_model.tar.gz"
    
    # Filter Mutect2 calls
    echo "Filtering Mutect2 variants for ${newname}..."
    gatk --java-options "${java_mem}" FilterMutectCalls \
        -R "${ref_genome}" \
        -V "${outputdir}/human_variants/${newname}_somatic_unfiltered.vcf.gz" \
        --ob-priors "${outputdir}/human_variants/${newname}_read_orientation_model.tar.gz" \
        -O "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz"

    # Generate variant statistics
    echo "Generating variant statistics..."
    bcftools stats "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz" > \
        "${outputdir}/reports/${newname}_variant_stats.txt"
    
    echo "✅ Mutect2 tumor-only variant calling completed for ${newname}"
    
    # Annotate variants with Funcotator (optional)
    echo "Annotating variants with Funcotator..."
    gatk --java-options "${java_mem}" Funcotator \
         --variant "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz" \
         --reference "${ref_genome}" \
         --output "${outputdir}/human_variants/${newname}_annotated.vcf.gz" \
         --output-file-format VCF \
         --data-sources-path ${funcotator_data} \
         --ref-version hg38
    
    # Generate variant statistics
    echo "Generating variant statistics..."
    bcftools stats "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz" > \
        "${outputdir}/reports/${newname}_variant_stats.txt"
    
    echo "✅ Mutect2 variant calling completed for ${newname}"

    ###############################################
    ###############################################
    ###############################################
    ###############################################

    # Extract reads directly mapped to HTLV-1
    echo "Extracting reads directly mapped to HTLV-1..."
    samtools view -@ ${nproc} -b "${bamq}_sorted.bam" "J02029.1" > \
    "${outputdir}/virus/${newname}_direct_viral.bam"

    # Sort and index
    samtools sort -@ ${nproc} -o "${outputdir}/virus/${newname}_direct_sorted.bam" \
    "${outputdir}/virus/direct/${newname}_direct_viral.bam"
    samtools index "${outputdir}/virus/${newname}_direct_sorted.bam"

    # Compare coverage statistics
    samtools coverage "${outputdir}/virus/direct/${newname}_direct_sorted.bam" > \
    "${outputdir}/reports/${newname}_direct_viral_coverage.txt"


   
    # Picard to mark duplicates
    ###$$$ Check after how to call this issue with picard
    echo "Marking duplicates in human reads..."
    java ${java_mem} -jar ${picard_cmd} MarkDuplicates  \
        -I "${outputdir}/virus/${newname}_direct_sorted.bam" \
        -O "${outputdir}/virus/${newname}_dedup.bam" \
        -M "${outputdir}/virus/${newname}_dedup_metrics.txt" \
        -CREATE_INDEX true \
        -REMOVE_DUPLICATES true \
        -VALIDATION_STRINGENCY LENIENT

  
    ######################################################
    # HTLV-1 VARIANT CALLING
    ######################################################
    echo "======== CALLING MUTATIONS IN HTLV-1 GENOME ========"

    # Call variants using GATK HaplotypeCaller tuned for stable virus
    echo "Calling HTLV-1 variants with GATK for ${newname}..."
    gatk --java-options "${java_mem}" HaplotypeCaller \
    -R "${htlv1_ref}" \
    -I "${outputdir}/virus/${newname}_dedup.bam" \
    -O "${outputdir}/virus_variants/${newname}_htlv1_variants_gatk.vcf.gz" \
    --min-pruning 1 \
    --min-dangling-branch-length 1 \
    --standard-min-confidence-threshold-for-calling 20.0 \
    --sample-ploidy 10 \
    --min-base-quality-score 25 \
    --annotation StrandBiasBySample \
    --pcr-indel-model NONE

    # Generate variant statistics
    echo "Generating HTLV-1 variant statistics..."
    bcftools stats "${outputdir}/virus_variants/${newname}_htlv1_variants_gatk.vcf.gz" > \
    "${outputdir}/reports/${newname}_htlv1_variants_stats.txt"

    echo "✅ HTLV-1 variant calling completed for ${newname}"
    
    # Clean up intermediate files
    rm -f "${outputdir}/virus/${newname}_direct_viral.bam" \
          "${outputdir}/virus/${newname}_direct_sorted.bam" 

    
done < "${project_list}"


# FIX MULTIQC:
echo "📊 Running MultiQC"
multiqc "${outputdir}/fastp" -o "${outputdir}/reports" -n "multiqc_report_postqc"
echo "✅ MultiQC report: ${outputdir}/multiqc_report_postqc.html"