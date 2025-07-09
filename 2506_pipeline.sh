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
java_mem="-Xmx8G -XX:MaxDirectMemorySize=16G"

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

if false
then 
### Notes
java --version
sudo update-alternatives --config java
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

fi

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
    if [ "${outputdir}/temp/${newname}_1.fq.gz" -ef "${sra}R1_001.fastq.gz" ] && [ "${outputdir}/temp/${newname}_2.fq.gz" -ef "${sra}R3_001.fastq.gz" ]; then
        echo "✅ Fastq files already copied for ${newname}. Skipping copy step."
    else
        echo "...Copying and renaming ${sra}R1_001.fastq.gz to ${newname}_1.fq.gz"
        rsync -av --inplace --no-whole-file "${sra}R1_001.fastq.gz" "${outputdir}/temp/${newname}_1.fq.gz"
        echo "...Copying and renaming ${sra}R3_001.fastq.gz to ${newname}_2.fq.gz"
        rsync -av --inplace --no-whole-file "${sra}R3_001.fastq.gz" "${outputdir}/temp/${newname}_2.fq.gz"
        echo "${newname}" >> "${outputdir}/${date_suffix}_foranalysis.txt"
        echo "${sra},${newname}" >> "${outputdir}/${date_suffix}_ref.txt"
        echo "✅ Fastq_1 and Fastq_2 of ${newname} copied to ${outputdir}/temp/"
    fi
    

    # Quality control with fastp
    #if [" ${outputdir}/temp/${newname}_postqc_1.fq.gz" && "${outputdir}/temp/${newname}_postqc_2.fq.gz" -f ]; then
    if [ -f "${outputdir}/temp/${newname}_postqc_1.fq.gz" ] && [ -f "${outputdir}/temp/${newname}_postqc_2.fq.gz" ]; then
        echo "✅ Post-QC files already exist for ${newname}. Skipping fastp step."
        continue
    fi
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
        --thread 16

    echo "✅ Quality control completed for ${newname}"
	
    
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
    samtools sort -@ ${nproc} -o "${bamq}_sorted.bam" 
    echo "✅ Sorted BAM file created: ${bamq}_sorted.bam"

    #INDEXING
    echo "Indexing BAM file"
    samtools index "${bamq}_sorted.bam"
    echo "✅ Indexed BAM file created: ${bamq}_sorted.bam.bai"

    


    ###############################################
    ###############################################
    ###############################################
    
    echo "======== EXTRACTING HUMAN READS OF INTEREST ========"

    # Extract well-mapped human reads (properly paired, MAPQ>=30)
    echo "Extracting high-quality human reads..."
    # Create a temporary BED file with human chromosomes
    samtools idxstats "${bamq}_sorted.bam" | grep -v "J02029.1" | awk '{print $1"\t0\t"$2}' > "${outputdir}/human/${newname}_human_regions.bed"

    # Use the BED file to extract human reads
    samtools view -@ ${nproc} -b -f 2 -q 30 -M -L "${outputdir}/human/${newname}_human_regions.bed" \
    "${bamq}_sorted.bam" > "${outputdir}/human/${newname}_human.bam"

    # Sort and index
    #samtools sort -@ ${nproc} -o "${outputdir}/human/${newname}_final.bam" "${outputdir}/human/${newname}_human.bam"
    #samtools index "${outputdir}/human/${newname}_final.bam"
    #echo "✅ Human reads extracted and sorted: ${outputdir}/human/${newname}_final.bam"
    
    ## Creating BED file for coverage regions > 30
    samtools depth -a "${outputdir}/human/${newname}_human.bam" | \
    awk '$3 >= 30 {print $1"\t"$2-1"\t"$2"\t"$3}' | \
    bedtools merge -d 10 -c 4 -o mean > "${outputdir}/human/${newname}_targeted.bed"
    echo "✅ High coverage regions BED file created: ${outputdir}/human/${newname}_targeted.bed"


    # Mark duplicates
    echo "Marking duplicates in human reads..."
    gatk --java-options "${java_mem}" MarkDuplicatesSpark \
    -I "${outputdir}/human/${newname}_human.bam" \
    -O "${outputdir}/human/${newname}_final_dedup.bam" \
    -M "${outputdir}/reports/${newname}_dedup_metrics.txt" \
    -L "${outputdir}/human/${newname}_targeted.bed" \
    --create-output-bam-index true \
    --remove-sequencing-duplicates true \
    --spark-master local[4] \
    --tmp-dir "${outputdir}/temp/"
    

    #java ${java_mem} -jar ${picard_cmd} MarkDuplicates \
    #    -I "${outputdir}/human/${newname}_final.bam" \
    #    -O "${outputdir}/human/${newname}_final_dedup.bam" \
    #    -M "${outputdir}/human/${newname}_dedup_metrics.txt" \
    #    -CREATE_INDEX true \
    #    -REMOVE_DUPLICATES true \
    #    -VALIDATION_STRINGENCY LENIENT

    
    echo "✅ Human reads processed and duplicates marked: ${outputdir}/human/${newname}_final_dedup.bam"
    

    # Remove intermediate files
    echo "Cleaning up intermediate files..."
    
    rm -f "${outputdir}/human/${newname}_final.bam" \
          "${outputdir}/human/${newname}_final.bam.bai" \
          "${outputdir}/human/${newname}_human.bam" \
          "${outputdir}/human/${newname}_human_regions.bed"
    rm -f "${fastq1}" "${fastq2}"
    rm -f "${outputdir}/temp/${newname}_1.fq.gz" \
        "${outputdir}/temp/${newname}_2.fq.gz"
    
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
     --ref-version hg38 \
     --transcript-selection-mode BEST_EFFECT \
     --exclude-field COSMIC_overlapping_mutations 
    
    
    # Generate variant statistics
    echo "Generating variant statistics..."
    bcftools stats "${outputdir}/human_variants/${newname}_somatic_filtered.vcf.gz" > \
        "${outputdir}/reports/${newname}_variant_stats.txt"
    
    echo "✅ Mutect2 variant calling completed for ${newname}"

    
    rm "${outputdir}/human/${newname}_final_dedup.bam" "${outputdir}/human/${newname}_final_dedup.bai"  
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
    "${outputdir}/virus/${newname}_direct_viral.bam"
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