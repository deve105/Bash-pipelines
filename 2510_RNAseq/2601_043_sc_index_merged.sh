#!/bin/bash


set -euo pipefail

echo "🧬 Indexing Fasta and GTF in STAR"

ref_dir="/home/htlvatl/coreapps/2025_refgenomes/2601_001_di_htlv1_ebv_hg38"
hg="GRCh38.primary_assembly.genome"
htlv1="AF033817_1"
ebv="AJ507799_2"

cd "$ref_dir"

mkdir -p "${ref_dir}/2601_001_di_htlv1_ebv_hg38"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1: Indexing using STAR"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"


STAR --runMode genomeGenerate \
echo "        --genomeDir star_index_merged \\"
echo "        --genomeFastaFiles $merged_fasta \\"
echo "        --sjdbGTFfile $merged_gtf \\"
echo "        --runThreadN \$(nproc) \\"
echo "        --sjdbOverhang 99"