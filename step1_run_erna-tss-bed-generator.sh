#!/bin/bash

# Input and output file paths
input_file="./data/enhancer_rna_files/enhancer_quartiles.csv"
output_file="./data/enhancer_rna_files/enhancer_quartiles.bed"

# Create the output file
> "$output_file"

# Skip the header line and process each row
tail -n +2 "$input_file" | while IFS=',' read -r index chromosome position mean_expression quartile; do
    # Format the BED line
    # Fields: chromosome, start, end, name, score, strand, group
    echo -e "chr$chromosome\t$((position-1))\t$position\tenhancer_region\t$mean_expression\t.\t$quartile" >> "$output_file"
done

echo "Conversion complete. Output saved to $output_file"
