#!/bin/bash


set -euo pipefail

echo "🧬 Merging all fasta and GTF"

ref_dir="/home/htlvatl/coreapps/2025_refgenomes"
hg="GRCh38.primary_assembly.genome"
htlv1="AF033817_1"
ebv="AJ507799_2"
cmv="AY446894_2"

cd "$ref_dir"

mkdir -p "${ref_dir}/2601_001_di_htlv1_ebv_hg38"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1: Merging FASTA Files"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "🔗 Merging fasta GRCh38 with HTLV-1, EBV and CMV..."

cat "${hg}.fa" "${htlv1}.fa" "${ebv}.fa" "${cmv}.fa" > "${ref_dir}/2601_001_di_htlv1_ebv_hg38/2601_001_rf_htlv1_ebv_hg38.fa"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1: Merging GTF Files"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "🔗 Merging GTF GRCh38 with HTLV-1 and EBV..."


merged_gtf="${ref_dir}/2601_001_di_htlv1_ebv_hg38/2601_001_rf_htlv1_ebv_hg38.gtf"
hg2="gencode.v49.primary_assembly.annotation"

echo "🔗 Merging GTF files..."

# Create merged GTF with headers and proper formatting
{
    # Copy human GTF header
    grep "^#" "${hg2}.gtf"
    
    # Add custom comment about merge
    echo "## HTLV-1 annotation merged for viral-host analysis"
    echo "## EBV annotation merged for viral-host analysis"
    
    # Copy human annotations (skip headers)
    grep -v "^#" "${hg2}.gtf"
    
    # Copy HTLV-1 annotations (skip headers)
    grep -v "^#" "${htlv1}.gtf"

    # Copy EBV annotations (skip headers)
    grep -v "^#" "${ebv}.gtf"

    # Copy CMV annotations (skip headers)
    grep -v "^#" "${cmv}.gtf"
    
} > "${merged_gtf}"

if [ $? -eq 0 ]; then
    echo "✅ Merged GTF created: $merged_gtf"
    echo "   Size: $(ls -lh "${merged_gtf}" | awk '{print $5}')"
else
    echo "❌ Error merging GTF files"
    exit 1
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 4: Cleaning Up Intermediate Files"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "🗑️  Removing intermediate files..."

# Method 1: Using array (recommended)
files_to_remove=(
    "${ref_dir}/${htlv1}.fa"
    "${ref_dir}/${htlv1}.gtf"
    "${ref_dir}/${ebv}.fa"
    "${ref_dir}/${ebv}.gtf"
    "${ref_dir}/${hg}.fa"
    "${ref_dir}/${hg2}.gtf"
    "${ref_dir}/${cmv}.fa"
    "${ref_dir}/${cmv}.gtf"
)

for file in "${files_to_remove[@]}"; do
    if [ -f "$file" ]; then
        rm -f "$file"
        echo "   ✓ Removed: $(basename "$file")"
    fi
done

echo "✅ Cleanup completed"