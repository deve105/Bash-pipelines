#!/bin/bash
# filepath: /home/htlvatl/Documents/GitHub/Bash-pipelines/2602_ENAdownload/2602_ENAdownloader.sh

set -euo pipefail  # Strict mode

echo "=========================================="
echo "Script for ENA downloader and retry logic"
echo "by Daniel Enriquez Version V1.4 2026-03-09"
echo "=========================================="
echo ""

# References
# https://github.com/enasequence/ena-ftp-downloader/
# https://ega-archive.org/datasets/EGAD00001015669
# https://www.ebi.ac.uk/ena/browser/view/PRJEB47382

# ============================================================================
# CONFIGURATION
# ============================================================================

disc_tmp="/home/htlvatl/Documents/devapps/temp"
hdd_tmp="/media/htlvatl/mem2"
downloader_ena="/home/htlvatl/coreapps/ena-file-downloader-1.1.8/ena-file-downloader.jar"

# Download retry settings
max_retries=5
retry_wait=30

# ============================================================================
# COLOR CODES FOR OUTPUT
# ============================================================================

readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly NC='\033[0m'

# ============================================================================
# LOGGING FUNCTIONS
# ============================================================================

log_error() {
    printf "%b[✗ ERROR]%b %s\n" "$RED" "$NC" "$*" >&2
}

log_success() {
    printf "%b[✓ SUCCESS]%b %s\n" "$GREEN" "$NC" "$*"
}

log_warn() {
    printf "%b[⚠ WARNING]%b %s\n" "$YELLOW" "$NC" "$*"
}

log_info() {
    printf "%b[ℹ INFO]%b %s\n" "$BLUE" "$NC" "$*"
}

# ============================================================================
# VALIDATION - FILE ARGUMENTS
# ============================================================================

# Check if argument is provided
if [ $# -eq 0 ]; then
    log_error "Usage: $0 <sample_list.txt>"
    log_info "Example: $0 PRJNA1292147.txt"
    exit 1
fi

# Check if file exists
if [ ! -f "$1" ]; then
    log_error "File not found: '$1'"
    log_info "Please provide a valid file with SRA/ENA accessions"
    exit 1
fi

# Check if file is non-empty
if [ ! -s "$1" ]; then
    log_error "File is empty: '$1'"
    exit 1
fi

log_success "Input file validated: $1"
echo ""

# ============================================================================
# VALIDATION - DEPENDENCIES
# ============================================================================

log_info "Checking dependencies..."

# Check Java installation
if ! command -v java &>/dev/null; then
    log_error "Java not found!"
    log_info "Install with: sudo apt install openjdk-8-jdk"
    exit 1
fi

# Check ENA downloader JAR
if [ ! -f "$downloader_ena" ]; then
    log_error "ENA downloader JAR not found: $downloader_ena"
    log_info "Build with: cd /path/to/ena-ftp-downloader && ./gradlew build"
    exit 1
fi

log_success "All dependencies found"
echo ""

# ============================================================================
# SETUP - FOLDER STRUCTURE
# ============================================================================

filename=$(basename "$1")
filename="${filename%.*}"  # Remove extension
filename_path="${hdd_tmp}/${filename}"
filename_disc="${disc_tmp}/${filename}"

log_info "Project name: $filename"
echo ""

# Step 1: Create main project folder
log_info "Step 1: Creating main project folder..."
if [ -d "$filename_path" ]; then
    log_warn "Folder already exists: $filename_path"
else
    mkdir -p "${filename_path}"/{rawdata,bam,logs,reports}
    log_success "Created: $filename_path"
fi
echo ""

# Step 2: Create temporary project folder
log_info "Step 2: Creating temporary project folder..."
if [ -d "$filename_disc" ]; then
    log_warn "Folder already exists: $filename_disc"
else
    mkdir -p "${filename_disc}"/{rawdata,bam,logs,reports}
    log_success "Created: $filename_disc"
fi
echo ""

# ============================================================================
# INITIALIZE VARIABLES
# ============================================================================

total_samples=$(wc -l < "$1")
failed_samples=()
successful_samples=()
skipped_samples=0

log_file="${filename_disc}/download.log"
: > "$log_file"

# ============================================================================
# CHECK DISK SPACE
# ============================================================================

log_info "Checking disk space..."
available_gb=$(df "$hdd_tmp" | tail -1 | awk '{printf int($4/1024/1024)}')
log_info "Available space on $hdd_tmp: ${available_gb}GB"

if [ "$available_gb" -lt 100 ]; then
    log_warn "Warning: Less than 100GB free!"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_error "Aborted"
        exit 1
    fi
fi
echo ""

# ============================================================================
# PROCESS SRA IDS WITH DOWNLOAD AND RETRY
# ============================================================================

log_info "Step 3: Processing SRA/ENA accessions..."
echo ""

while IFS= read -r sra; do
    # Skip empty lines
    [[ -z "$sra" ]] && ((skipped_samples++)) && continue
    
    # Skip comments
    [[ "$sra" =~ ^# ]] && ((skipped_samples++)) && continue
    
    # Trim whitespace
    sra=$(echo "$sra" | xargs)
    
    # Check if already downloaded (avoid re-downloading)
    output_dir="${filename_path}/rawdata"
    if [ -d "$output_dir" ] && find "$output_dir" -name "*${sra}*" -type f | grep -q .; then
        log_success "  ✓ Already exists: $sra (skipping)"
        successful_samples+=("$sra")
        continue
    fi
    
    # Download with retry logic
    success=false
    last_error=""
    
    for attempt in $(seq 1 "$max_retries"); do
        log_info "  Attempt $attempt/$max_retries: $sra"
        
        # Capture error output for debugging
        if java -jar "${downloader_ena}" \
            --location="${output_dir}" \
            --accessions="${sra}" \
            --format=READS_FASTQ \
            --protocol=FTP \
            --asperaLocation=null >> "$log_file" 2>&1; then
            
            success=true
            log_success "  ✓ Downloaded: $sra"
            successful_samples+=("$sra")
            break
            
        else
            last_error="$?"
            
            if [ $attempt -lt "$max_retries" ]; then
                # Exponential backoff: wait longer with each retry
                wait_time=$((retry_wait * attempt))
                log_warn "  Attempt $attempt failed (exit code: $last_error)"
                log_warn "  Retrying in ${wait_time}s..."
                sleep "$wait_time"
            fi
        fi
    done
    
    if [ "$success" = false ]; then
        log_error "  ✗ Failed: $sra (exit code: $last_error)"
        failed_samples+=("$sra")
        echo "[FAILED] $sra (error code: $last_error)" >> "$log_file"
    fi
    
done < "$1"
