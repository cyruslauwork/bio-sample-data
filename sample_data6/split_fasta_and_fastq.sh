#!/bin/bash

# Input FASTA and FASTQ files
input_fasta="src/1000Mbp_MP_merged/1000Mbp_MP_merged_m8_SortedBySStart_MSA_degap"
input_fastq="src/1000Mbp_MP_merged/1000Mbp_MP_merged_m8_SortedBySStart_MSA_degap.fastq"

# Temporary file to store shuffled sequence headers
temp_headers="temp_headers.txt"

# Shuffle sequence headers (for pairing purposes)
grep '^>' "$input_fasta" | shuf > "$temp_headers"

# Function to extract sequences from FASTA based on headers
extract_fasta_sequences() {
  awk 'BEGIN{RS=">"; ORS=""} NR==FNR{headers[$1]; next} $1 in headers {print ">"$0}' "$1" "$2"
}

# Function to extract sequences from FASTQ based on headers
extract_fastq_sequences() {
  awk 'NR==FNR{headers[$0]; next} $1 in headers {print; getline; print; getline; print; getline; print}' "$1" "$2"
}

# Create files of approximate sizes
sizes=(10 25 50 100) # Sizes in MB
for size in "${sizes[@]}"; do
  # Calculate approximate number of sequences needed
  seq_count=$(($(wc -l < "$temp_headers") * size * 1048576 / 143000000))
  
  # Select random headers and extract sequences for both FASTA and FASTQ
  head -n "$seq_count" "$temp_headers" > "selected_headers.txt"
  extract_fasta_sequences "selected_headers.txt" "$input_fasta" > "output_${size}MB.fasta"
  sed 's/^>//' "selected_headers.txt" > "selected_headers_no_gt.txt"
  extract_fastq_sequences "selected_headers_no_gt.txt" "$input_fastq" > "output_${size}MB.fastq"
done

# Cleanup
rm "$temp_headers" "selected_headers.txt" "selected_headers_no_gt.txt"