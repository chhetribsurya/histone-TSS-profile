#!/bin/bash

# This script runs the gencode_tss_bed_generator.py script to generate TSS BED files 
# from Gencode annotations and expression group files

# Default settings
OUTPUT_DIR="results/gencode_tss_bed"
#EXPR_DIR="data/expression_groups"  # Directory with expression groups files
EXPR_DIR="./results/expression_quantiles_with0"  # Directory with expression groups files
GENOME="hg19"
GENCODE_VERSION="37"  # Use numeric version (without "v" prefix)

FILE_PATTERN="ENS"    # Pattern in filenames (e.g., "TCGA_GTEX_genes_ENS_0.txt")
GENE_PREFIX="ENSG"    # Prefix for gene IDs
TSS_SELECTION="canonical"

# Create output directory
mkdir -p $OUTPUT_DIR

# Basic usage without RefSeq matching
#python gencode_tss_bed_generator.py \
#  --expr_dir $EXPR_DIR \
#  --output_dir $OUTPUT_DIR \
#  --output_bed "gencode_tss_with_groups.bed" \
#  --genome $GENOME \
#  --gencode_version $GENCODE_VERSION \
#  --tss_selection $TSS_SELECTION \
#  --file_pattern $FILE_PATTERN \
#  --gene_prefix $GENE_PREFIX

# To enable RefSeq matching, uncomment and run this instead:
#python gencode_tss_bed_generator.py \
python gencode_tss_bed_generator_v2.py \
   --expr_dir $EXPR_DIR \
   --output_dir $OUTPUT_DIR \
   --output_bed "gencode_tss_with_groups.bed" \
   --genome $GENOME \
   --gencode_version $GENCODE_VERSION \
   --tss_selection $TSS_SELECTION \
   --match_refseq \
   --file_pattern $FILE_PATTERN \
   --gene_prefix $GENE_PREFIX

echo "Gencode TSS BED files have been generated in $OUTPUT_DIR"
