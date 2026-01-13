#!/bin/bash


set -euo pipefail

echo "🧬 Indexing Fasta in STAR"

ref_dir="/home/htlvatl/coreapps/2025_refgenomes/2601_001_di_htlv1_ebv_hg38"
merged_fasta="${ref_dir}/2601_001_rf_htlv1_ebv_hg38.fa"
merged_gtf="${ref_dir}/2601_001_rf_htlv1_ebv_hg38.gtf"
index_dir="${ref_dir}/2601_001_rf_star_index_merged"
overhang=99 
#SAindex=min(14, log2(GenomeLength)/2 - 1)
SAindex=14

cd "$ref_dir"

# Determine number of threads (use nproc for auto-detection)
#num_threads=${1:-$(nproc)}
num_threads=10
# Calculate RAM needed (STAR recommends ~30GB per genome copy)
# For safety with large genomes, use 60GB

genome_size_limit=60000000000

# ============================================================================
# VALIDATION
# ============================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1: Validation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Check if FASTA file exists
if [ ! -f "$merged_fasta" ]; then
    echo "❌ Error: FASTA file not found: $merged_fasta"
    exit 1
fi
echo "✅ FASTA file found: $(basename "$merged_fasta") ($(ls -lh "$merged_fasta" | awk '{print $5}'))"

# Check if GTF file exists
if [ ! -f "$merged_gtf" ]; then
    echo "❌ Error: GTF file not found: $merged_gtf"
    exit 1
fi
echo "✅ GTF file found: $(basename "$merged_gtf") ($(ls -lh "$merged_gtf" | awk '{print $5}'))"

# Check if FASTA has index
if [ ! -f "${merged_fasta}.fai" ]; then
    echo "⚠️  FASTA index not found. Creating..."
    if command -v samtools &>/dev/null; then
        samtools faidx "$merged_fasta"
        echo "✅ FASTA index created"
    else
        echo "⚠️  samtools not found, continuing without index"
    fi
fi

# Check if STAR is installed
if ! command -v STAR &>/dev/null; then
    echo "❌ Error: STAR is not installed or not in PATH"
    echo "   Install with: mamba install -c bioconda star"
    exit 1
fi
echo "✅ STAR found: $(STAR --version)"

# Check available resources
echo ""
echo "📊 System Resources:"
echo "   Available CPU threads: $(nproc)"
echo "   Using threads: $num_threads"
if [ -f /proc/meminfo ]; then
    total_mem=$(awk '/MemTotal/ {print int($2/1024/1024) "GB"}' /proc/meminfo)
    echo "   Total memory: $total_mem"
fi

# ============================================================================
# CLEANUP OLD INDEX (if exists)
# ============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 2: Preparing Index Directory"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -d "$index_dir" ]; then
    echo "⚠️  Index directory already exists: $index_dir"
    read -p "Remove and recreate? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$index_dir"
        echo "🗑️  Old index directory removed"
    else
        echo "⚠️  Using existing directory (may cause issues)"
    fi
fi

mkdir -p "$index_dir"
echo "✅ Index directory ready: $index_dir"

# Change to working directory
cd "$ref_dir"



# ============================================================================
# STEP 3: GENERATE STAR INDEX
# ============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 3: Generating STAR Index"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "⏱️  This may take 15-30 minutes depending on genome size..."
echo ""

# Generate STAR index
if STAR --runMode genomeGenerate \
    --genomeDir "$index_dir" \
    --genomeFastaFiles "$merged_fasta" \
    --sjdbGTFfile "$merged_gtf" \
    --runThreadN "$num_threads" \
    --sjdbOverhang "$overhang" \
    --limitGenomeGenerateRAM "$genome_size_limit" \
    --outTmpDir "${index_dir}_tmp" \
    --outFileNamePrefix "${index_dir}/" \
    --genomeSAindexNbases "${SAindex}/"; then
    
    echo ""
    echo "✅ STAR index generation completed successfully!"
    
    # Clean up temporary directory
    if [ -d "${index_dir}_tmp" ]; then
        rm -rf "${index_dir}_tmp"
        echo "🧹 Temporary files cleaned up"
    fi
else
    echo ""
    echo "❌ Error: STAR index generation failed!"
    echo "   Check the STAR log files in: $index_dir"
    exit 1
fi

# ============================================================================
# STEP 4: VALIDATION
# ============================================================================


echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 4: Validating Index"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Check for required index files
required_files=(
    "SA"
    "SAindex"
    "Genome"
    "genomeParameters.txt"
    "chrLength.txt"
    "chrNameLength.txt"
    "chrName.txt"
    "chrStart.txt"
)

echo "🔍 Checking index files..."
all_files_present=true

for file in "${required_files[@]}"; do
    if [ -f "${index_dir}/${file}" ]; then
        size=$(ls -lh "${index_dir}/${file}" | awk '{print $5}')
        echo "   ✅ $file ($size)"
    else
        echo "   ❌ $file (MISSING)"
        all_files_present=false
    fi
done

if [ "$all_files_present" = true ]; then
    echo "✅ All required index files present"
else
    echo "❌ Some index files are missing!"
    exit 1
fi


# ============================================================================
# STEP 5: SUMMARY
# ============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ STAR INDEX GENERATION COMPLETED"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo ""
echo "📊 Index Statistics:"
echo "   Genome FASTA: $(basename "$merged_fasta")"
echo "   Annotation GTF: $(basename "$merged_gtf")"
echo "   Index directory: $index_dir"
echo "   Index size: $(du -sh "$index_dir" | awk '{print $1}')"
echo ""
echo "📁 Index files:"
ls -lh "$index_dir" | tail -n +2 | awk '{printf "   - %-30s (%5s)\n", $9, $5}'
echo ""
echo "🎯 Ready for alignment! Use this index with STAR mapping:"
echo "   STAR --genomeDir $index_dir \\"
echo "        --readFilesIn sample_R1.fastq sample_R2.fastq \\"
echo "        --outFileNamePrefix sample_ \\"
echo "        --runThreadN $num_threads"
echo ""