#!/bin/bash

set -euo pipefail

echo "HLA Pipeline Remote Deployment Script"
echo "Version 1.0 - $(date)"

# Configuration - EDIT THESE VALUES
REMOTE_HOST=""  # e.g., "user@server.example.com"
REMOTE_BASE="/home/user/bioinformatics"  # Remote working directory

# Local paths
LOCAL_SCRIPTS="/home/htlvatl/Documents/GitHub/Bash-pipelines"

# Function to deploy pipeline
deploy_pipeline() {
    echo "=== Deploying HLA Pipeline to Remote Host ==="
    
    # Create remote directory structure
    ssh "${REMOTE_HOST}" "
        mkdir -p ${REMOTE_BASE}/{scripts,data,output,logs,references}
        mkdir -p ${REMOTE_BASE}/Peru_IRID/{fastq,bam,reports,HLA-LA}
    "
    
    # Copy scripts
    echo "Copying pipeline scripts..."
    scp "${LOCAL_SCRIPTS}/2503_peruremapping.sh" "${REMOTE_HOST}:${REMOTE_BASE}/scripts/"
    scp "${LOCAL_SCRIPTS}/2506_pipeline.sh" "${REMOTE_HOST}:${REMOTE_BASE}/scripts/"
    scp "${LOCAL_SCRIPTS}/install_bwa_mem2.sh" "${REMOTE_HOST}:${REMOTE_BASE}/scripts/"
    
    # Create environment setup script
    echo "Creating remote environment setup..."
    ssh "${REMOTE_HOST}" "cat > ${REMOTE_BASE}/setup_environment.sh << 'EOF'
#!/bin/bash

echo 'Setting up bioinformatics environment...'

# Update package manager
sudo apt update

# Install essential packages
sudo apt install -y build-essential git wget curl

# Install bioinformatics tools via conda/mamba
if ! command -v conda &> /dev/null; then
    echo 'Installing Miniforge (includes conda and mamba)...'
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh -b -p \$HOME/miniforge3
    export PATH=\$HOME/miniforge3/bin:\$PATH
    echo 'export PATH=\$HOME/miniforge3/bin:\$PATH' >> ~/.bashrc
fi

# Create bioinformatics environment
mamba create -n biotools -c bioconda -c conda-forge \\
    bwa-mem2 samtools fastp multiqc python=3.9 -y

echo 'Environment setup complete!'
echo 'Activate with: conda activate biotools'
EOF"
    
    chmod +x "${REMOTE_HOST}:${REMOTE_BASE}/setup_environment.sh"
}

# Function to execute pipeline remotely
execute_hla_pipeline() {
    local INPUT_DIR="$1"
    
    echo "=== Executing HLA Pipeline Remotely ==="
    
    # Update script paths for remote execution
    ssh "${REMOTE_HOST}" "
        cd ${REMOTE_BASE}
        
        # Activate conda environment
        source ~/miniforge3/etc/profile.d/conda.sh
        conda activate biotools
        
        # Update script variables for remote paths
        sed -i 's|/home/labatl/devprojects|${REMOTE_BASE}|g' scripts/2503_peruremapping.sh
        sed -i 's|/home/htlvatl/Documents/devapps|${REMOTE_BASE}|g' scripts/2503_peruremapping.sh
        
        # Start pipeline in background
        nohup bash scripts/2503_peruremapping.sh > logs/hla_pipeline_\$(date +%Y%m%d_%H%M%S).log 2>&1 &
        
        echo 'Pipeline started!'
        echo 'Monitor with: tail -f ${REMOTE_BASE}/logs/hla_pipeline_*.log'
    "
}

# Function to monitor pipeline
monitor_pipeline() {
    echo "=== Monitoring Remote Pipeline ==="
    
    ssh "${REMOTE_HOST}" "
        cd ${REMOTE_BASE}
        
        echo 'Running processes:'
        ps aux | grep -E '(bwa-mem2|samtools|fastp|HLA-LA)' | grep -v grep
        
        echo ''
        echo 'Latest log:'
        ls -t logs/*.log 2>/dev/null | head -1 | xargs tail -20 || echo 'No logs found'
        
        echo ''
        echo 'Disk usage:'
        du -sh * 2>/dev/null
    "
}

# Function to download results
download_results() {
    echo "=== Downloading Results ==="
    
    mkdir -p results_$(date +%Y%m%d)
    
    rsync -avz --progress "${REMOTE_HOST}:${REMOTE_BASE}/output/" "results_$(date +%Y%m%d)/"
    rsync -avz --progress "${REMOTE_HOST}:${REMOTE_BASE}/logs/" "results_$(date +%Y%m%d)/logs/"
    
    echo "Results downloaded to: results_$(date +%Y%m%d)/"
}

# Main execution
case "${1:-}" in
    "check")
        ./check_remote_system.sh "${REMOTE_HOST}"
        ;;
    "deploy")
        deploy_pipeline
        ;;
    "setup")
        ssh "${REMOTE_HOST}" "bash ${REMOTE_BASE}/setup_environment.sh"
        ;;
    "execute")
        execute_hla_pipeline "${2:-}"
        ;;
    "monitor")
        monitor_pipeline
        ;;
    "download")
        download_results
        ;;
    *)
        echo "Usage: $0 {check|deploy|setup|execute|monitor|download}"
        echo ""
        echo "Commands:"
        echo "  check    - Check remote system specifications"
        echo "  deploy   - Deploy pipeline scripts to remote host"
        echo "  setup    - Setup bioinformatics environment"
        echo "  execute  - Execute HLA pipeline"
        echo "  monitor  - Monitor running pipeline"
        echo "  download - Download results"
        echo ""
        echo "First, edit this script to set REMOTE_HOST variable"
        ;;
esac
