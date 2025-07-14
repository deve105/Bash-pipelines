#!/bin/bash

awk '
BEGIN {
    OFS = "\t"
    FS = "\t"
    first_file = 1
    in_header = 1
}

# Process each file
FNR == 1 {
    # Extract base filename for sample ID
    split(FILENAME, parts, "/")
    sample_id = parts[length(parts)]
    sub(/\..*$/, "", sample_id)  # Remove extension
    sub(/(_final|_maf)?$/, "", sample_id)  # Remove common suffixes
    in_header = 1
}

# Process header/metadata lines
in_header && /^#/ {
    if (first_file) {
        # Find Tumor_Sample_Barcode column in header
        if ($0 !~ /^##/ && $0 ~ /Tumor_Sample_Barcode/) {
            for (i=1; i<=NF; i++) {
                if ($i == "Tumor_Sample_Barcode") {
                    tumor_col = i
                    break
                }
            }
            if (!tumor_col) tumor_col = 16  # Fallback column
        }
        print $0
    }
    next
}

# Mark end of header
in_header && !/^#/ {
    in_header = 0
    first_file = 0
}

# Process data rows
!in_header {
    # Assign sample ID to Tumor_Sample_Barcode column
    if (tumor_col <= NF) {
        $tumor_col = sample_id
    } else {
        # Add column if missing
        $0 = $0 OFS sample_id
    }
    print $0
}
' *.maf > merged.maf

# Verify output
if [ -s merged.maf ]; then
    echo "Successfully merged MAF files"
    echo "First few lines:"
    head -n 5 merged.maf
    echo -e "\nSample distribution:"
    awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="Tumor_Sample_Barcode"){col=i;break}} 
                 NR>1 && col && !/^#/ {print $col}' merged.maf | sort | uniq -c
else
    echo "Error: Failed to create merged MAF file" >&2
    exit 1
fi