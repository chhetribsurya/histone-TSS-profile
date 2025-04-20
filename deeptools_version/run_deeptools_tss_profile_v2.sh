#!/usr/bin/bash

#SBATCH --job-name=plot_tss_profile
#SBATCH --output=./results/tss_profile_%A_%a.out
#SBATCH --error=./results/tss_profile_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=6
#SBATCH --partition=normal
# Note: No --array parameter here, it will be provided at submission time

# Base output directory
BASE_OUTPUT_DIR="./results"
mkdir -p ${BASE_OUTPUT_DIR}

# Path to the script
SCRIPT_PATH="./deeptools_tss_profile_v8.py"

# Define WIG file to use
WIG_FILE="/data/baca/users/sc1238/datasets/scripts/for_alexis/data/GSM7578579_H1.3_WT.wig.gz"

# Create an array of all gene TSV files
GENE_FILES=($(ls /data/baca/users/sc1238/datasets/scripts/for_alexis/quantile_results_new_normalized/equal_bins/RNASeq*.tsv))

# Get the current file from the array based on SLURM_ARRAY_TASK_ID
CURRENT_GENE_FILE="${GENE_FILES[$SLURM_ARRAY_TASK_ID]}"

# Extract basename without extension for naming the output directory
BASENAME=$(basename "$CURRENT_GENE_FILE" .tsv)
OUTPUT_DIR="${BASE_OUTPUT_DIR}/${BASENAME}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Default parameter values
CONDA_ENV="genometools_env"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

# Print job info
echo "Starting job at $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID, Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing file: $CURRENT_GENE_FILE"
echo "Output directory: $OUTPUT_DIR"

# Run the script with specified WIG file and gene list
python $SCRIPT_PATH \
  --wig ${WIG_FILE} \
  --genes ${CURRENT_GENE_FILE} \
  --gtf_version 37 \
  --genome hg19 \
  --groups 1,2,3,4,5,6,7,8,9,10,All \
  --output_prefix $(basename "$WIG_FILE" .wig.gz) \
  --output_dir ${OUTPUT_DIR}

echo "Job completed at $(date)"
