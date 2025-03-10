#!/bin/bash

#set #-euo pipefail # Strict mode

echo "Script to concatenate txt files from HLA-LA output"
echo "Version V1.3 2025-03-10"

maindir="/home/labatl/devprojects/Peru_IRID"
output_file="${maindir}/hla_la_concatenated_output.txt"
skipped_log="${maindir}/hla_la_skipped_log.txt"

# Clear the output file and skipped files log
> "$output_file"
> "$skipped_log"

# Check if input file exists
if [[ ! -f "${maindir}/Peru_IRID.txt" ]]; then
    echo "Error: Input file 'Peru_IRID.txt' not found in $maindir" >&2
    exit 1
fi

# Counter to track which sample we're processing
sample_counter=0

while IFS= read -r sample; do
    # Increment the sample counter
    ((sample_counter++))
    
    filex="${maindir}/HLA-LA/${sample}/hla/R1_bestguess_G.txt"
    
    # Check if file exists
    if [[ ! -f "$filex" ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping '$filex' (file not found)" >> "$skipped_log"
        continue
    fi
    
    # Count the number of lines in the file
    line_count=$(wc -l < "$filex")
    
    # Skip the file if it is empty or has only one line
    if [[ "$line_count" -le 1 ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping '$filex' (empty or only one row)" >> "$skipped_log"
        continue
    fi
    
    # If this is the first sample in the loop, add the header row
    if [[ $sample_counter -eq 1 ]]; then
        echo "Adding header from first sample: $sample"
        head -n 1 "$filex" | awk '{print $0 "\tSample"}' >> "$output_file"
    fi
    
    # Append data to output file, adding sample name as last column
    tail -n +2 "$filex" | awk -v sample="${sample}" '{print $0 "\t" sample}' >> "$output_file"
done < "${maindir}/Peru_IRID.txt"

echo "Processing complete. Results in $output_file"
echo "Skipped files logged in $skipped_log"