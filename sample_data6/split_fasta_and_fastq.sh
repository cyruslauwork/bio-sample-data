#!/bin/bash

# Input FASTA file
input_fasta="src/1000Mbp_MP_merged/1000Mbp_MP_merged_m8_SortedBySStart_MSA_degap"

# Temporary file to store shuffled sequence headers
temp_headers="temp_headers.txt"

# Shuffle sequence headers
grep '^>' "$input_fasta" | shuf > "$temp_headers"

# Function to extract sequences based on headers
extract_sequences() {
  awk 'BEGIN{RS=">"; ORS=""} NR==FNR{headers[$1]; next} $1 in headers {print ">"$0}' "$1" "$2"
}

# Create files of approximate sizes
sizes=(10 25 50 100) # Sizes in MB
for size in "${sizes[@]}"; do
  # Calculate approximate number of sequences needed
  seq_count=$(($(wc -l < "$temp_headers") * size * 1048576 / 143000000))
  
  # Select random headers and extract sequences
  head -n "$seq_count" "$temp_headers" > "selected_headers.txt"
  extract_sequences "selected_headers.txt" "$input_fasta" > "output_${size}MB.fasta"
done

# Cleanup
rm "$temp_headers" "selected_headers.txt"