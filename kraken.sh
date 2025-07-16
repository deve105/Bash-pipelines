#!/bin/bash


# Deltaretrovirus Genome Fetcher and Kraken2 DB Preprocessor

# 1. Set up directories
mkdir -p ./sequences

# 2. Define core deltaretrovirus TaxIDs (update as needed)
cat > taxids.txt <<EOF
11901
11902
11903
11904
11905
11906
11907
11908
11909
11926
11927
11928
1454033
153136
1907191
194440
194441
194443
2772270
28332
3128424
318275
318279
33747
33748
3428211
3428212
3428213
3428214
347956
36368
36369
37814
39015
39016
39101
402036
402042
402043
402044
402046
406769
47688
481147
EOF

# 3. Fetch and process each genome
while read taxid; do
    [[ "$taxid" =~ ^# ]] && continue  # Skip comments
    
    echo "Processing TaxID: $taxid"
    
    # A. Download genome assembly
    datasets download genome taxon $taxid --filename ./sequences/${taxid}_tmp.zip
    
    # B. Extract and rename sequences
    unzip -j ./sequences/${taxid}_tmp.zip -d sequences/${taxid}
    rm ./sequences/${taxid}_tmp.zip
    
    # Process each FASTA file
    for f in sequences/${taxid}/*.fna; do
        # Extract accession (first word after >)
        acc=$(grep "^>" "$f" | head -1 | cut -d" " -f1 | tr -d ">")
        
        # Create standardized name: taxid_accession.fna
        new_name="sequences/${taxid}_${acc}.fna"
        
        # Process headers to Kraken2 format: >accession|kraken:taxid|taxid
        sed "s/^>.*$/>${acc}|kraken:taxid|${taxid}/" "$f" > "$new_name"
        
        # Remove original file
        rm "$f"
    done
done < taxids.txt

# 4. Create sequence to TaxID mapping (required for Kraken2)
find sequences -name "*.fna" | while read f; do
    acc=$(basename "$f" .fna | cut -d_ -f2)
    taxid=$(basename "$f" .fna | cut -d_ -f1)
    grep "^>" "$f" | awk -v taxid="$taxid" '{print $1 "\t" taxid}' | tr -d ">"
done > seqid2taxid.map
if false; then
# 5. Build Kraken2 database
kraken2-build --add-to-library sequences/*.fna --db deltaretro_db
kraken2-build --build --db deltaretro_db \
    --threads 8 \
    --kmer-len 25 \
    --minimizer-len 20 \
    --minimizer-spaces 5
fi
echo "Database built in: kraken2_db/deltaretro_db"