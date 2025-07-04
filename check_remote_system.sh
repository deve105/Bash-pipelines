#!/bin/bash

# Quick remote system checker
# Usage: ./check_remote_system.sh user@hostname

if [ $# -eq 0 ]; then
    echo "Usage: $0 user@hostname"
    echo "Example: $0 labuser@192.168.1.100"
    exit 1
fi

REMOTE_HOST="$1"

echo "Connecting to ${REMOTE_HOST}..."
echo "======================================"

ssh "${REMOTE_HOST}" '
echo "🖥️  SYSTEM INFORMATION"
echo "Hostname: $(hostname)"
echo "OS: $(cat /etc/os-release | grep PRETTY_NAME | cut -d'"'"'\"'"'"' -f2 2>/dev/null || uname -s)"
echo "Kernel: $(uname -r)"
echo "Uptime: $(uptime | cut -d'"'"','"'"' -f1 | cut -d'"'"' '"'"' -f4-)"
echo ""

echo "🔧 CPU INFORMATION"
echo "CPU Model: $(lscpu | grep '"'"'Model name'"'"' | cut -d'"'"':'"'"' -f2 | sed '"'"'s/^[ \t]*//'"'"')"
echo "CPU Cores: $(nproc) cores"
echo "CPU Threads: $(lscpu | grep '"'"'^CPU(s):'"'"' | awk '"'"'{print $2}'"'"')"
echo "Architecture: $(uname -m)"
echo ""

echo "💾 MEMORY INFORMATION"
free -h | grep -E '"'"'Mem|Swap'"'"'
echo ""

echo "💿 DISK SPACE"
df -h | head -5
echo ""

echo "🧬 BIOINFORMATICS SOFTWARE"
echo "Checking available software..."
for tool in bwa-mem2 bwa samtools fastp multiqc python3 conda mamba; do
    if command -v $tool >/dev/null 2>&1; then
        VERSION=$(${tool} --version 2>&1 | head -1 | cut -d'"'"' '"'"' -f1-2 2>/dev/null || echo "installed")
        echo "✅ $tool: $VERSION"
    else
        echo "❌ $tool: not found"
    fi
done
echo ""

echo "🐍 PYTHON ENVIRONMENT"
if command -v python3 >/dev/null 2>&1; then
    echo "Python: $(python3 --version)"
    echo "Pip packages related to bioinformatics:"
    pip3 list 2>/dev/null | grep -E '"'"'(bio|genomics|seq|alignment)'"'"' | head -5 || echo "No bioinformatics packages found"
fi
echo ""

echo "📊 SYSTEM LOAD"
echo "Load average: $(uptime | awk -F'"'"'load average:'"'"' '"'"'{print $2}'"'"')"
echo "Active processes: $(ps aux | wc -l)"
'
