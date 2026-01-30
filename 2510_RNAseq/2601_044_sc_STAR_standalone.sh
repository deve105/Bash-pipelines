#!/bin/bash

set -euo pipefail #Stric mode

echo "Script for STAR mapping with overhang 149" 
echo "by Daniel Enriquez Version V1.2 2026-01-07"


platform="ILLUMINA"
numthreads=16
overhang=149

fastq1=$1
fastq2=$2
sra=$3
ref_dir="/home/htlvatl/coreapps/2025_refgenomes/2601_001_di_htlv1_hg38_ADJH"
merged_fasta="/home/htlvatl/coreapps/2025_refgenomes/2601_001_di_htlv1_hg38_ADJH/2601_001_rf_htlv1_hg38_ADJH.fa"
merged_gtf="/home/htlvatl/coreapps/2025_refgenomes/2601_001_di_htlv1_hg38_ADJH/2601_001_rf_htlv1_hg38_ADJH.gtf"


# STAR index specific to overhang
index_dir="${ref_dir}/2601_001_rf_star_index_merged_h149_ADJH"

STAR --twopassMode Basic \
                --runThreadN 12 \
                --genomeDir "${index_dir}" \
                --readFilesIn "${fastq1}" "${fastq2}" \
                --readFilesCommand zcat \
                --sjdbOverhang "${overhang}" \
                --outSAMtype BAM SortedByCoordinate \
                --sjdbGTFfile "${merged_gtf}"  \
                --twopass1readsN -1 \
                --quantMode TranscriptomeSAM GeneCounts \
                --outSAMattrRGline ID:$sra SM:$sra PL:$platform \
                --outFileNamePrefix "${sra}_2" \
                --outFilterMultimapScoreRange 1 \
                --outFilterMultimapNmax 20 \
                --outFilterMismatchNmax 10 \
                --alignIntronMax 500000 \
                --alignMatesGapMax 1000000 \
                --sjdbScore 2 \
                --alignSJDBoverhangMin 1 \
  	            --outFilterMatchNminOverLread 0.33 \
                --outFilterScoreMinOverLread 0.33 \
                --outSAMattributes NH HI NM MD AS XS nM \
                --outSAMunmapped Within \
                --chimSegmentMin 12 \
                --chimJunctionOverhangMin 8 \
                --chimOutJunctionFormat 1 \
                --alignSJstitchMismatchNmax 5 -1 5 5 \
                --chimMultimapScoreRange 3 \
                --chimScoreJunctionNonGTAG -4 \
                --chimMultimapNmax 20 \
                --chimNonchimScoreDropMin 10 \
                --peOverlapNbasesMin 12 \
                --peOverlapMMp 0.1 \
                --alignInsertionFlush Right \
                --alignSplicedMateMapLminOverLmate 0 \
                --alignSplicedMateMapLmin 30