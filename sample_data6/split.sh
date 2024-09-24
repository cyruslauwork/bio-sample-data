#!/bin/bash

sam_file="output.sam"
fasta_file="igv_fasta.fa"
bam_split_prefix="results/split_bam/split_"
fasta_split_prefix="results/split_fa/split_"
sequences_per_file=10000  # Adjust this number as needed

# Create output directories
mkdir -p results/split_bam results/split_fa

# Function to split FASTA file, create ID mapping, and index split FASTAs
split_fasta() {
    awk -v seqs_per_file="$sequences_per_file" -v prefix="$fasta_split_prefix" '
    BEGIN {part=1; count=0; file=prefix sprintf("%03d", part) ".fa"}
    /^>/ {
        if (count >= seqs_per_file) {
            close(file);
            part++;
            file=prefix sprintf("%03d", part) ".fa";
            count=0;
        }
        if (count > 0) {
            print "" > file;
        }
        count++;
        id = $0;
        sub(/^>/, "", id);
        print id "\t" sprintf("%03d", part) > "fasta_id_mapping.txt";
    }
    {print > file}' "$fasta_file"

    # Index each split FASTA file
    for fa_file in ${fasta_split_prefix}*.fa; do
        samtools faidx "$fa_file"
    done
}

# Function to split SAM file based on FASTA splits and create indexed BAM files
split_sam_to_bam() {
    # Extract header from SAM file
    grep '^@' "$sam_file" > sam_header.sam

    # Process SAM file and split
    grep -v '^@' "$sam_file" | \
    awk -v prefix="$bam_split_prefix" '
    BEGIN {
        while ((getline < "fasta_id_mapping.txt") > 0) {
            split($0, a, "\t");
            map[a[1]] = a[2];
        }
        close("fasta_id_mapping.txt");
    }
    {
        if ($3 in map) {
            outfile = prefix map[$3] ".sam";
            print > outfile;
        }
    }'

    # Convert split SAM files to BAM, add headers, sort, and index
    for split_sam_file in ${bam_split_prefix}*.sam; do
        bam_file="${split_sam_file%.sam}.bam"
        sorted_bam_file="${bam_file%.bam}.sorted.bam"
        cat sam_header.sam "$split_sam_file" | samtools view -bS - | samtools sort -o "$sorted_bam_file"
        samtools index "$sorted_bam_file"
        rm "$split_sam_file"
        mv "$sorted_bam_file" "$bam_file"
        mv "${sorted_bam_file}.bai" "${bam_file}.bai"
    done

    # Clean up
    rm sam_header.sam
}

# Split FASTA file, create ID mapping, and index split FASTAs
split_fasta

# Split SAM file based on FASTA splits, convert to BAM, and index split BAMs
split_sam_to_bam

# Clean up temporary files
rm fasta_id_mapping.txt

echo "Splitting and indexing complete. FASTA and BAM files are split into pairs with corresponding .fai and .bai files."