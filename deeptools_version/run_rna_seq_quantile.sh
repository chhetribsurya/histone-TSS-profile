#!/usr/bin/bash

INPUT_FILE="/data/baca/users/sc1238/datasets/scripts/for_alexis/data/GSE236538_RNASeq_H1Xsh_raw_read_counts_EnsemblID.txt"
OUTPUT_DIR="./quantile_results_new_normalized_V4"

# Default parameter values
CONDA_ENV="py_env"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

# Print job info
echo "Starting job at $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Execute the script
# python rna_seq_quantile_v5.py --input $INPUT_FILE --output_dir $OUTPUT_DIR --equal_bins
# python rna_seq_quantile_v7.py --input $INPUT_FILE --output_dir $OUTPUT_DIR --normalize --equal_bins

# TPM Normalization Only
# Download GENCODE data if needed
# Normalize the RNA-seq counts to TPM values
# Skip all quantile calculations
# python rna_seq_quantile_v8.py --input $INPUT_FILE --output_dir $OUTPUT_DIR --normalize_matrix_only

# Full Analysis (Default)
# Perform TPM normalization
# Calculate quantiles for each condition
# Create separate quantile files for each condition
python rna_seq_quantile_v8.py --input $INPUT_FILE --output_dir $OUTPUT_DIR --normalize --equal_bins

# Print completion message
echo "Job completed at $(date)"

