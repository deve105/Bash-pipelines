#!/bin/bash

# Loop through all files that match the pattern *.bam
for file in *.bam; do
    # Check if the file exists (in case there are no matches)
    if [[ -e "$file" ]]; then
        # Extract the new name by removing the prefix (e.g., "17_")
        newname=$(echo "$file" | sed 's/^[0-9]*_//')
        
        # Rename the .bam file
        mv "$file" "$newname"
        echo "Renamed '$file' to '$newname'"

        # Check if the corresponding .bam.bai file exists and rename it
        if [[ -e "${file}.bai" ]]; then
            mv "${file}.bai" "${newname}.bai"
            echo "Renamed '${file}.bai' to '${newname}.bai'"
        fi
    else
        echo "No files found matching the pattern *.bam"
    fi
done

# Generate a list of the renamed .bam files (without the .bam extension)
ls | grep -E "\.bam$" | sed 's/\.bam$//' > Peru_IRID.txt
echo "Generated list of renamed files in Peru_IRID.txt"