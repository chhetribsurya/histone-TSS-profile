#!/usr/bin/bash

# Default input data directory and file
DATA_DIR='/data/baca/users/ss2976/GeneExp_Project/M2GEP_Shared_Repo/Data/Other_Data_Sources/'
DEFAULT_FILE_NAME='All_Genes_TCGA_GTex_RSEM_TPM_converted.rds'
DEFAULT_FILE_PATH="${DATA_DIR}/${DEFAULT_FILE_NAME}"

# Create output directory if it doesn't exist
OUTPUT_DIR="results/expression_quantiles_with0"
mkdir -p $OUTPUT_DIR

# Set parameters
NUM_GROUPS=10
PREFIX="TCGA_GTEX_genes"
ID_COLUMN="gene_id_no_version"

# Default parameter values
CONDA_ENV="r_bio_env"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

# To include non-expressed genes as Group 0, uncomment these lines:
INCLUDE_NON_EXPRESSED="--include-non-expressed"
NON_EXPRESSED_THRESHOLD=0.01
NON_EXPRESSED_THRESHOLD=0.05

# Build command
CMD="Rscript expression_quantiles_v2.R \
  --input \"$DEFAULT_FILE_PATH\" \
  --output-dir \"$OUTPUT_DIR\" \
  --num-groups $NUM_GROUPS \
  --prefix \"$PREFIX\" \
  --id-column \"$ID_COLUMN\""

# Add non-expressed options if enabled
if [ ! -z ${INCLUDE_NON_EXPRESSED+x} ]; then
  CMD="$CMD $INCLUDE_NON_EXPRESSED --non-expressed-threshold $NON_EXPRESSED_THRESHOLD"
fi

# Execute the command
echo "Running: $CMD"
eval $CMD

echo "Expression quantile groups have been saved to $OUTPUT_DIR"
