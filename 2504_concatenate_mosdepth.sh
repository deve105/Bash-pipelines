#!/bin/bash

#set #-euo pipefail # Strict mode

echo "Script to concatenate txt files from HLA-LA output"
echo "Version V1.3 2025-03-10"


for file in *_v2.regions.bed.gz; do
    filename=$(basename "$file" _v2.regions.bed.gz)  # Removes .gz extension
    zcat "$file" | awk -v fname="$filename" '{print $0 "\t" fname}' >> ../2504_regions_combined_output.tsv
done