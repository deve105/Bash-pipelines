#!/bin/bash
# filepath: /home/htlvatl/Documents/GitHub/Bash-pipelines/2510_RNAseq/create_gtf_from_ncbi.sh

set -euo pipefail

echo "🧬 Creating GTF file from NCBI data"
echo "Supports both GenBank and RefSeq annotations"

# Function to download NCBI data and convert to GTF
create_gtf_from_ncbi() {
    local accession="$1"
    local organism_name="${2:-unknown}"
    local output_dir="${3:-ncbi_annotations}"
    
    echo "📥 Processing NCBI accession: $accession"
    
    # Create output directory
    mkdir -p "$output_dir"
    cd "$output_dir"
    
    # Download GenBank file
    echo "📥 Downloading GenBank file for $accession..."
    wget -O "${accession}.gb" \
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=${accession}&db=nuccore&report=genbank&retmode=text"
    
    # Download FASTA sequence
    echo "📥 Downloading FASTA sequence for $accession..."
    wget -O "${accession}.fasta" \
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=${accession}&db=nuccore&report=fasta&retmode=text"
    
    # Convert GenBank to GTF using Python script
    echo "🔄 Converting GenBank to GTF format..."
    
    python3 << EOF
import re
import sys
from datetime import datetime

def parse_genbank_to_gtf(gb_file, accession, organism):
    """Convert GenBank file to GTF format"""
    
    with open(gb_file, 'r') as f:
        content = f.read()
    
    # Extract basic information
    definition_match = re.search(r'DEFINITION\s+(.+?)(?=\n[A-Z])', content, re.DOTALL)
    definition = definition_match.group(1).replace('\n', ' ').strip() if definition_match else "Unknown"
    
    # Extract sequence length
    length_match = re.search(r'LOCUS\s+\S+\s+(\d+)\s+bp', content)
    seq_length = int(length_match.group(1)) if length_match else 0
    
    gtf_file = f"{accession}.gtf"
    
    with open(gtf_file, 'w') as out:
        # Write GTF header
        out.write('##gff-version 2\n')
        out.write(f'##date: {datetime.now().strftime("%Y-%m-%d")}\n')
        out.write(f'##genome-build: {accession}\n')
        out.write(f'##organism: {organism}\n')
        out.write(f'##sequence-length: {seq_length}\n')
        out.write(f'##definition: {definition}\n')
        out.write('#\n')
        
        # Find features section
        features_match = re.search(r'FEATURES\s+Location/Qualifiers\n(.*?)(?=ORIGIN|CONTIG|//)', content, re.DOTALL)
        if not features_match:
            print("⚠️  No features found in GenBank file")
            return gtf_file
        
        features_text = features_match.group(1)
        
        # Parse individual features
        feature_blocks = re.split(r'\n     (?=[a-zA-Z])', features_text)
        
        gene_counter = 0
        transcript_counter = 0
        
        for block in feature_blocks:
            if not block.strip():
                continue
                
            lines = block.split('\n')
            if not lines:
                continue
                
            # Parse feature line
            feature_line = lines[0].strip()
            feature_match = re.match(r'(\S+)\s+(.+)', feature_line)
            if not feature_match:
                continue
                
            feature_type = feature_match.group(1)
            location = feature_match.group(2)
            
            # Parse location
            if '..' in location:
                # Handle various location formats
                location_clean = re.sub(r'[<>]', '', location)
                location_clean = re.sub(r'complement\((.+)\)', r'\1', location_clean)
                
                if '..' in location_clean:
                    try:
                        start, end = location_clean.split('..')
                        start = int(re.sub(r'[^\d]', '', start))
                        end = int(re.sub(r'[^\d]', '', end))
                        
                        # Determine strand
                        strand = '-' if 'complement' in location else '+'
                        
                    except (ValueError, IndexError):
                        continue
                else:
                    continue
            else:
                continue
            
            # Parse qualifiers
            qualifiers = {}
            for line in lines[1:]:
                line = line.strip()
                if not line.startswith('/'):
                    continue
                    
                qual_match = re.match(r'/(\w+)=?"?([^"]*)"?', line)
                if qual_match:
                    key, value = qual_match.groups()
                    qualifiers[key] = value.strip('"')
            
            # Extract common attributes
            gene_name = qualifiers.get('gene', qualifiers.get('locus_tag', f'gene_{gene_counter + 1}'))
            product = qualifiers.get('product', qualifiers.get('note', 'unknown'))
            protein_id = qualifiers.get('protein_id', '')
            
            # Generate GTF entries based on feature type
            if feature_type == 'gene':
                gene_counter += 1
                gene_id = gene_name
                
                # Gene entry
                out.write(f'{accession}\tNCBI\tgene\t{start}\t{end}\t.\t{strand}\t.\t')
                out.write(f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "protein_coding"; ')
                out.write(f'gene_description "{product}";\n')
                
            elif feature_type == 'CDS':
                transcript_counter += 1
                gene_id = gene_name
                transcript_id = f"{gene_name}_001"
                
                # Gene entry (if not already created)
                out.write(f'{accession}\tNCBI\tgene\t{start}\t{end}\t.\t{strand}\t.\t')
                out.write(f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "protein_coding"; ')
                out.write(f'gene_description "{product}";\n')
                
                # Transcript entry
                out.write(f'{accession}\tNCBI\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t')
                out.write(f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_name "{gene_name}"; ')
                out.write(f'gene_type "protein_coding"; transcript_type "protein_coding";\n')
                
                # Exon entry
                out.write(f'{accession}\tNCBI\texon\t{start}\t{end}\t.\t{strand}\t.\t')
                out.write(f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "1"; ')
                out.write(f'gene_name "{gene_name}"; gene_type "protein_coding"; transcript_type "protein_coding";\n')
                
                # CDS entry
                out.write(f'{accession}\tNCBI\tCDS\t{start}\t{end}\t.\t{strand}\t0\t')
                out.write(f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "1"; ')
                out.write(f'gene_name "{gene_name}"; gene_type "protein_coding"; transcript_type "protein_coding"; ')
                if protein_id:
                    out.write(f'protein_id "{protein_id}"; ')
                out.write(f'product "{product}";\n')
                
            elif feature_type in ['rRNA', 'tRNA', 'ncRNA', 'misc_RNA']:
                gene_counter += 1
                transcript_counter += 1
                gene_id = gene_name
                transcript_id = f"{gene_name}_001"
                
                # Determine gene type
                if feature_type == 'rRNA':
                    gene_type = 'rRNA'
                elif feature_type == 'tRNA':
                    gene_type = 'tRNA'
                else:
                    gene_type = 'ncRNA'
                
                # Gene entry
                out.write(f'{accession}\tNCBI\tgene\t{start}\t{end}\t.\t{strand}\t.\t')
                out.write(f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}"; ')
                out.write(f'gene_description "{product}";\n')
                
                # Transcript entry
                out.write(f'{accession}\tNCBI\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t')
                out.write(f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_name "{gene_name}"; ')
                out.write(f'gene_type "{gene_type}"; transcript_type "{gene_type}";\n')
                
                # Exon entry
                out.write(f'{accession}\tNCBI\texon\t{start}\t{end}\t.\t{strand}\t.\t')
                out.write(f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "1"; ')
                out.write(f'gene_name "{gene_name}"; gene_type "{gene_type}"; transcript_type "{gene_type}";\n')
    
    return gtf_file

# Run the conversion
gtf_file = parse_genbank_to_gtf("${accession}.gb", "${accession}", "${organism_name}")
print(f"✅ GTF file created: {gtf_file}")
EOF
    
    echo "✅ GTF conversion completed"
    cd ..
}

# Function to create HTLV-1 GTF specifically
create_htlv1_gtf_from_ncbi() {
    echo "🦠 Creating HTLV-1 GTF from NCBI..."
    
    # HTLV-1 complete provirus
    create_gtf_from_ncbi "NC_001436.1" "Human T-lymphotropic virus 1" "htlv1_ncbi"
    
    # Alternative accessions for HTLV-1
    # create_gtf_from_ncbi "J02029.1" "Human T-lymphotropic virus 1" "htlv1_ncbi_alt"
    
    echo "✅ HTLV-1 GTF creation completed"
}

# Function to validate and enhance GTF
enhance_gtf() {
    local gtf_file="$1"
    local enhanced_gtf="${gtf_file%.gtf}_enhanced.gtf"
    
    echo "🔧 Enhancing GTF file: $gtf_file"
    
    # Add additional annotations and fix format issues
    awk '
    BEGIN {
        FS = OFS = "\t"
    }
    
    /^#/ {
        print $0
        next
    }
    
    NF >= 9 {
        # Ensure proper GTF format
        if ($3 == "gene") {
            # Add gene biotype if missing
            if ($9 !~ /gene_biotype/) {
                $9 = $9 " gene_biotype \"protein_coding\";"
            }
        }
        
        if ($3 == "transcript") {
            # Add transcript biotype if missing
            if ($9 !~ /transcript_biotype/) {
                $9 = $9 " transcript_biotype \"protein_coding\";"
            }
        }
        
        print $0
    }
    ' "$gtf_file" > "$enhanced_gtf"
    
    echo "✅ Enhanced GTF created: $enhanced_gtf"
}

# Function to merge with human GTF
merge_with_human() {
    local viral_gtf="$1"
    local human_gtf="$2"
    local merged_gtf="merged_human_viral.gtf"
    
    echo "🔗 Merging viral GTF with human annotation..."
    
    if [ ! -f "$human_gtf" ]; then
        echo "⚠️  Human GTF not found: $human_gtf"
        echo "Downloading human GTF..."
        wget -O gencode.v49.primary_assembly.annotation.gtf.gz \
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz"
        gunzip gencode.v49.primary_assembly.annotation.gtf.gz
        human_gtf="gencode.v49.primary_assembly.annotation.gtf"
    fi
    
    # Merge GTF files
    {
        # Human GTF headers
        grep "^#" "$human_gtf"
        echo "# Viral annotation from NCBI merged"
        echo "#"
        
        # Human annotations
        grep -v "^#" "$human_gtf"
        
        # Viral annotations
        grep -v "^#" "$viral_gtf"
        
    } > "$merged_gtf"
    
    echo "✅ Merged GTF created: $merged_gtf"
    echo "📊 Statistics:"
    echo "   Total entries: $(wc -l < "$merged_gtf")"
    echo "   Human genes: $(grep -c "chr[0-9XYM]" "$merged_gtf" || echo "0")"
    echo "   Viral genes: $(grep -c "NC_001436.1\|J02029.1" "$merged_gtf" || echo "0")"
}

# Main function
main() {
    case "${1:-}" in
        --htlv1)
            create_htlv1_gtf_from_ncbi
            ;;
        --merge)
            if [ $# -lt 3 ]; then
                echo "Usage: $0 --merge <viral.gtf> <human.gtf>"
                exit 1
            fi
            merge_with_human "$2" "$3"
            ;;
        --enhance)
            if [ $# -lt 2 ]; then
                echo "Usage: $0 --enhance <file.gtf>"
                exit 1
            fi
            enhance_gtf "$2"
            ;;
        *)
            if [ $# -lt 1 ]; then
                echo "Usage: $0 <NCBI_accession> [organism_name] [output_dir]"
                echo "       $0 --htlv1"
                echo "       $0 --merge <viral.gtf> <human.gtf>"
                echo "       $0 --enhance <file.gtf>"
                echo ""
                echo "Examples:"
                echo "  $0 NC_001436.1 'Human T-lymphotropic virus 1'"
                echo "  $0 NC_045512.2 'SARS-CoV-2'"
                echo "  $0 --htlv1"
                exit 1
            fi
            create_gtf_from_ncbi "$1" "${2:-unknown}" "${3:-ncbi_annotations}"
            ;;
    esac
}

# Run main function
main "$@"