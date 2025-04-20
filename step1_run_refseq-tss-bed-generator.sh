#!/bin/bash

# Run RefSeq TSS BED Generator
# This script downloads RefSeq gene annotations and creates TSS BED files 
# with expression groups from the provided gene lists

# Define variables
OUTPUT_DIR="refseq_tss_results"
EXPR_DIR="./data/T47D_ExpressionGroups_RNASeq"
GENOME="hg19"
#GENOME="hg38"

# Default parameter values
CONDA_ENV="genometools_env"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

# Create output directory
mkdir -p $OUTPUT_DIR

# Print job info
echo "Starting job at $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Run the generator script
echo "Starting RefSeq TSS BED Generator..."
#python refseq-tss-bed-generator_v2.py \
python refseq-tss-bed-generator.py \
  --expr_dir $EXPR_DIR \
  --output_dir $OUTPUT_DIR \
  --genome $GENOME

echo "Done!"
