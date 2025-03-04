#!/bin/bash

set -euo pipefail #Stric mode

project_directory="/home/labatl/devprojects/Peru_IRID/optitype"

while IFS= read -r sra; do
    echo "HLA typing for ${sra}"
    fastq_1="${sra}_sorted_chr6_1.fastq"
    fastq_2="${sra}_sorted_chr6_2.fastq"
    sudo docker run \
        -v "${project_directory}:/data/" \
        -t fred2/optitype \
        -d \
        -i "/data/${fastq_1}" "/data/${fastq_2}" \
        -o "/data/output/" \
        -p "${sra}_"
done < "$1"



