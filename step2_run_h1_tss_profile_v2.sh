#!/bin/bash

# H1 histone variant distribution analysis script
# Analyzes H1.0 distribution around TSS grouped by expression levels

# Define paths
#OUTPUT_DIR="./resultsFINAL/results_H1.3.tcgaGtexGenes"
OUTPUT_DIR="./resultsFINAL/results_H1.0.eRNA-quartiles"
OUTPUT_DIR="./resultsFINAL/results_H1.2.eRNA-quartiles"
OUTPUT_DIR="./resultsFINAL/results_H1.3.eRNA-quartiles"

OUTPUT_FILE="${OUTPUT_DIR}/h1_tss_profile.png"
LOG_FILE="${OUTPUT_DIR}/analysis_log.txt"

TSS_FILE="./refseq_tss_results/refseq_tss_with_groups.bed"
TSS_FILE="./scripts_stable_genes_version/results/gencode_tss_bed/matched_refseq_tss_with_groups.bed"
TSS_FILE="./scripts_stable_genes_version/results/gencode_tss_bed/gencode_tss_with_groups.bed"
TSS_FILE="./data/enhancer_rna_files/combined_enhancer_eRNA_tss_with_groups.bed"
TSS_FILE="./data/enhancer_rna_files/enhancer_quartiles.bed"

#WIG_FILE="./data//GSM5076926_H1.0_WT.wig"
#WIG_FILE="./data/GSM5076927_H1.2_WT.wig"
WIG_FILE="./data/GSM7578579_H1.3_WT.wig"
#WIG_FILE="./data/GSM5076928_H1.4_WT.wig"
#WIG_FILE="./data/GSM5076929_H1.5_WT.wig"
#WIG_FILE="./data/GSM5076930_H1X_WT.wig"

# Parameters
WINDOW_SIZE=3000
BIN_SIZE=10
PROCESSES=6                 # Number of CPU cores to use
#EXPRESSION_GROUPS="0,1,2,3,4,5,6,7,8,9,10"  # Exclude group 0 (no expression)
EXPRESSION_GROUPS="1,2,3,4"  # Exclude group 0 (no expression)
CONVERT_TO_BIGWIG=true      # Convert WIG to BigWig for faster processing
CHROM_SIZES="hg19.chrom.sizes"  # Required for WIG to BigWig conversion
#SAMPLE_SIZE=10000           # Set to empty to process all TSS positions

# Default parameter values
CONDA_ENV="genometools_env"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Print start time
echo "Starting H1 histone variant analysis at $(date)" | tee "$LOG_FILE"
echo "TSS file: $TSS_FILE" | tee -a "$LOG_FILE"
echo "WIG file: $WIG_FILE" | tee -a "$LOG_FILE"
echo "Output will be saved to: $OUTPUT_FILE" | tee -a "$LOG_FILE"
echo "-------------------------------------------" | tee -a "$LOG_FILE"

# Check if input files exist
if [ ! -f "$TSS_FILE" ]; then
    echo "Error: TSS file $TSS_FILE not found!" | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -f "$WIG_FILE" ]; then
    echo "Error: WIG file $WIG_FILE not found!" | tee -a "$LOG_FILE"
    exit 1
fi

# Build command
#--abs
#--log2 \
#--format gff \
#--input $GFF_FILE \

echo "Running H1 distribution analysis..."
#CMD="python h1_tss_profile_v3.py \
CMD="python h1_tss_profile.py \
  --tss $TSS_FILE \
  --input $WIG_FILE \
  --format wig \
  --output $OUTPUT_FILE \
  --window $WINDOW_SIZE \
  --bin $BIN_SIZE \
  --processes $PROCESSES"

# Add optional parameters
if [ -n "$EXPRESSION_GROUPS" ]; then
    CMD="$CMD --groups $EXPRESSION_GROUPS"
fi

if [ -n "$SAMPLE_SIZE" ]; then
    CMD="$CMD --sample $SAMPLE_SIZE"
fi

if [ "$CONVERT_TO_BIGWIG" = true ]; then
    CMD="$CMD --convert --chrom-sizes $CHROM_SIZES"
fi

# Run the analysis script
echo "Command: $CMD"
$CMD
echo "Task completed..."


