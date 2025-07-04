#!/bin/bash

set -euo pipefail

echo "Remote Pipeline Execution Script"
echo "Version 1.0 - $(date)"

# Configuration
REMOTE_HOST=""  # e.g., "user@192.168.1.100" or "user@server.domain.com"
REMOTE_WORKDIR=""  # e.g., "/home/user/bioinformatics"
LOCAL_PIPELINE_DIR="/home/htlvatl/Documents/GitHub/Bash-pipelines"

# Function to check remote system specs
check_remote_specs() {
    echo "=== Checking Remote System Specifications ==="
    ssh "${REMOTE_HOST}" "
        echo 'Hostname: \$(hostname)'
        echo 'OS: \$(cat /etc/os-release | grep PRETTY_NAME | cut -d'\"' -f2)'
        echo ''
        echo '=== CPU Information ==='
        lscpu | grep -E '^CPU\(s\)|^Model name|^Thread|^Core'
        echo ''
        echo '=== Memory Information ==='
        free -h
        echo ''
        echo '=== Disk Space ==='
        df -h | head -5
        echo ''
        echo '=== Available Software ==='
        which bwa-mem2 || echo 'bwa-mem2: not found'
        which samtools || echo 'samtools: not found'
        which fastp || echo 'fastp: not found'
        which multiqc || echo 'multiqc: not found'
    "
}

# Function to setup remote environment
setup_remote_environment() {
    echo "=== Setting up Remote Environment ==="
    
    # Create remote working directory
    ssh "${REMOTE_HOST}" "mkdir -p ${REMOTE_WORKDIR}/{scripts,data,output,logs}"
    
    # Copy pipeline scripts
    echo "Copying pipeline scripts..."
    scp "${LOCAL_PIPELINE_DIR}"/*.sh "${REMOTE_HOST}:${REMOTE_WORKDIR}/scripts/"
    
    # Copy any reference files (if small enough)
    echo "Note: Large reference files should be downloaded directly on remote host"
}

# Function to execute pipeline remotely
execute_remote_pipeline() {
    local SCRIPT_NAME="$1"
    local LOG_FILE="pipeline_$(date +%Y%m%d_%H%M%S).log"
    
    echo "=== Executing ${SCRIPT_NAME} on Remote Host ==="
    
    # Execute with logging and keep session alive
    ssh "${REMOTE_HOST}" "
        cd ${REMOTE_WORKDIR}
        nohup bash scripts/${SCRIPT_NAME} > logs/${LOG_FILE} 2>&1 &
        echo 'Pipeline started with PID: \$!'
        echo 'Log file: ${REMOTE_WORKDIR}/logs/${LOG_FILE}'
        echo 'Monitor with: tail -f ${REMOTE_WORKDIR}/logs/${LOG_FILE}'
    "
}

# Function to monitor remote pipeline
monitor_remote_pipeline() {
    echo "=== Monitoring Remote Pipeline ==="
    ssh "${REMOTE_HOST}" "
        cd ${REMOTE_WORKDIR}
        echo 'Active processes:'
        ps aux | grep -E '(bwa-mem2|samtools|fastp|multiqc)' | grep -v grep
        echo ''
        echo 'Latest log entries:'
        ls -t logs/*.log | head -1 | xargs tail -20
    "
}

# Function to sync results back
sync_results() {
    echo "=== Syncing Results Back ==="
    mkdir -p "${LOCAL_PIPELINE_DIR}/remote_results"
    rsync -avz "${REMOTE_HOST}:${REMOTE_WORKDIR}/output/" "${LOCAL_PIPELINE_DIR}/remote_results/"
    rsync -avz "${REMOTE_HOST}:${REMOTE_WORKDIR}/logs/" "${LOCAL_PIPELINE_DIR}/remote_results/logs/"
}

# Main menu
main_menu() {
    echo ""
    echo "=== Remote Pipeline Execution Menu ==="
    echo "1. Check remote system specs"
    echo "2. Setup remote environment"
    echo "3. Execute HLA mapping pipeline"
    echo "4. Execute virus detection pipeline"
    echo "5. Monitor running pipeline"
    echo "6. Sync results back"
    echo "7. Interactive SSH session"
    echo "8. Exit"
    echo ""
    read -p "Choose option (1-8): " choice
    
    case $choice in
        1) check_remote_specs ;;
        2) setup_remote_environment ;;
        3) execute_remote_pipeline "2503_peruremapping.sh" ;;
        4) execute_remote_pipeline "2506_pipeline.sh" ;;
        5) monitor_remote_pipeline ;;
        6) sync_results ;;
        7) ssh "${REMOTE_HOST}" ;;
        8) exit 0 ;;
        *) echo "Invalid option" ;;
    esac
}

# Check if configuration is set
if [[ -z "$REMOTE_HOST" ]]; then
    echo "Please configure REMOTE_HOST and REMOTE_WORKDIR variables first"
    echo "Edit this script and set:"
    echo "REMOTE_HOST=\"user@hostname.or.ip\""
    echo "REMOTE_WORKDIR=\"/path/to/remote/workdir\""
    exit 1
fi

# Run main menu
while true; do
    main_menu
    echo ""
    read -p "Press Enter to continue..."
done
