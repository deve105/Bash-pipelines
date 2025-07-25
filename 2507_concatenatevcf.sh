first=1
for file in *_breakend.vcf; do
    filename=$(basename "$file" _breakend.vcf)
    if [ $first -eq 1 ]; then
        # Include the first row from the first file
        cat "$file" | grep -v "^##" | awk -v fname="$filename" '{print $0 "\t" fname}' >> ../2507_concatenatedtab_vcf.tsv
        first=0
    else
        # Exclude the first row from subsequent files
        cat "$file" | grep -v "^##" | tail -n +2 | awk -v fname="$filename" '{print $0 "\t" fname}' >> ../2507_concatenatedtab_vcf.tsv
    fi
done