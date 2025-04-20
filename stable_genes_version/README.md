# Stable Genes Analysis Pipeline

This directory contains scripts for analyzing gene expression data and generating TSS (Transcription Start Site) BED files. The pipeline consists of two main components:

1. Expression Quantile Analysis
2. TSS BED File Generation

## 1. Expression Quantile Analysis

### Overview
The expression quantile analysis scripts process gene expression data to group genes into quantiles based on their expression levels. This helps in identifying genes with similar expression patterns across different conditions.

### Scripts

#### `expression_quantiles_v2.R`
Main R script for calculating expression quantiles. It processes gene expression data and groups genes into specified quantiles.

**Key Features:**
- Processes gene expression data from RDS files
- Groups genes into user-defined number of quantiles
- Option to include non-expressed genes as a separate group
- Handles both healthy (GTEX) and cancer (TCGA) samples
- Generates detailed output files for each quantile group

**Usage:**
```bash
Rscript expression_quantiles_v2.R \
  --input <input_file> \
  --output-dir <output_directory> \
  --num-groups <number_of_quantiles> \
  --prefix <output_prefix> \
  --id-column <gene_id_column> \
  [--include-non-expressed] \
  [--non-expressed-threshold <threshold>]
```

**Parameters:**
- `--input`: Path to input RDS file containing gene expression data
- `--output-dir`: Directory to save output files
- `--num-groups`: Number of quantile groups to create
- `--prefix`: Prefix for output files
- `--id-column`: Column name containing gene IDs
- `--include-non-expressed`: Flag to include non-expressed genes
- `--non-expressed-threshold`: Threshold for considering genes as non-expressed

#### `run_expression_quantiles_v2.sh`
Wrapper script for running the expression quantile analysis.

**Default Settings:**
- Input data directory: `/data/baca/users/ss2976/GeneExp_Project/M2GEP_Shared_Repo/Data/Other_Data_Sources/`
- Default input file: `All_Genes_TCGA_GTex_RSEM_TPM_converted.rds`
- Output directory: `results/expression_quantiles_with0`
- Number of groups: 10
- Prefix: "TCGA_GTEX_genes"
- ID column: "gene_id_no_version"

## 2. TSS BED File Generation

### Overview
The TSS BED file generation scripts create BED files containing transcription start site information for genes, optionally grouped by expression levels.

### Scripts

#### `gencode_tss_bed_generator_v2.py`
Python script for generating TSS BED files from Gencode annotations and expression group files.

**Key Features:**
- Processes Gencode annotation files
- Supports multiple genome versions (hg19, hg38)
- Option to match with RefSeq annotations
- Supports different TSS selection methods
- Generates BED files with expression group information

**Usage:**
```bash
python gencode_tss_bed_generator_v2.py \
  --expr_dir <expression_directory> \
  --output_dir <output_directory> \
  --output_bed <output_bed_file> \
  --genome <genome_version> \
  --gencode_version <gencode_version> \
  --tss_selection <selection_method> \
  [--match_refseq] \
  --file_pattern <pattern> \
  --gene_prefix <prefix>
```

**Parameters:**
- `--expr_dir`: Directory containing expression group files
- `--output_dir`: Directory to save output BED files
- `--output_bed`: Name of output BED file
- `--genome`: Genome version (e.g., hg19, hg38)
- `--gencode_version`: Gencode version number
- `--tss_selection`: TSS selection method (e.g., "canonical")
- `--match_refseq`: Flag to enable RefSeq matching
- `--file_pattern`: Pattern in filenames
- `--gene_prefix`: Prefix for gene IDs

#### `run_gencode_tss_bed_generator.sh`
Wrapper script for running the TSS BED file generation.

**Default Settings:**
- Output directory: `results/gencode_tss_bed`
- Expression directory: `./results/expression_quantiles_with0`
- Genome version: hg19
- Gencode version: 37
- File pattern: "ENS"
- Gene prefix: "ENSG"
- TSS selection: "canonical"

## Workflow

1. Run expression quantile analysis:
   ```bash
   ./run_expression_quantiles_v2.sh
   ```

2. Generate TSS BED files:
   ```bash
   ./run_gencode_tss_bed_generator.sh
   ```

## Dependencies

- R (with packages: optparse, dplyr, matrixStats)
- Python
- Conda environment: r_bio_env

## Output Files

### Expression Quantile Analysis
- `TCGA_GTEX_genes_ENS_<group>.txt`: Files containing genes in each quantile group
- `TCGA_GTEX_genes_ENS_0.txt`: File containing non-expressed genes (if enabled)

### TSS BED Generation
- `gencode_tss_with_groups.bed`: BED file containing TSS information with expression groups

## Notes

- The pipeline is designed to work with TCGA and GTEX gene expression data
- Expression values are expected to be in TPM (Transcripts Per Million) format
- The scripts support both hg19 and hg38 genome versions
- RefSeq matching is optional and can be enabled when needed

## Contact

For questions, improvement or issues, please contact: chhetribsurya@gmail.com
