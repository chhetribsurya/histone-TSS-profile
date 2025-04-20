# DeepTools-based TSS Profile Analysis

This directory contains scripts for analyzing transcription start site (TSS) profiles using deepTools, an alternative implementation to the main pipeline. This version provides additional features and optimizations for large-scale analysis.

## Overview

The deepTools version implements three main analysis steps:

1. **RNA-seq Quantile Analysis** (`rna_seq_quantile_v8.py`)
   - Processes RNA-seq data
   - Performs TPM normalization
   - Calculates expression quantiles
   - Groups genes by expression levels

2. **TSS Profile Generation** (`deeptools_tss_profile_v8.py`)
   - Uses deepTools for profile generation
   - Processes WIG/BigWig files
   - Creates TSS profiles for different expression groups
   - Supports parallel processing

3. **Stable Gene Analysis** (`stable_genes.R`)
   - Identifies stable genes across conditions
   - Performs statistical analysis
   - Generates gene lists for downstream analysis

## Scripts and Usage

### 1. RNA-seq Quantile Analysis

```bash
# Basic usage
./run_rna_seq_quantile.sh

# Direct script usage
python rna_seq_quantile_v8.py \
  --input /path/to/RNAseq_data.txt \
  --output_dir ./quantile_results \
  --normalize \
  --equal_bins
```

#### Parameters

| Parameter | Description | Default | Required | Example |
|-----------|-------------|---------|----------|---------|
| `--input` | Input RNA-seq file | None | Yes | `GSE236538_RNASeq_H1Xsh_raw_read_counts_EnsemblID.txt` |
| `--output_dir` | Output directory | None | Yes | `./quantile_results` |
| `--normalize` | Perform TPM normalization | `False` | No | `True` |
| `--equal_bins` | Use equal-sized bins | `False` | No | `True` |
| `--normalize_matrix_only` | Only normalize, skip quantiles | `False` | No | `True` |
| `--n_bins` | Number of quantile bins | `10` | No | `5` |

### 2. TSS Profile Generation

```bash
# Basic usage
./run_deeptools_tss_profile_v2.sh

# Direct script usage
python deeptools_tss_profile_v8.py \
  --wig /path/to/signal.wig \
  --genes /path/to/gene_list.tsv \
  --gtf_version 37 \
  --genome hg19 \
  --groups 1,2,3,4,5 \
  --output_prefix H1.3 \
  --output_dir ./results
```

#### Parameters

| Parameter | Description | Default | Required | Example |
|-----------|-------------|---------|----------|---------|
| `--wig` | Input WIG/BigWig file | None | Yes | `GSM7578579_H1.3_WT.wig.gz` |
| `--genes` | Gene list file | None | Yes | `RNASeq_quantile_1.tsv` |
| `--gtf_version` | GTF version | `37` | No | `37` or `38` |
| `--genome` | Genome version | `hg19` | No | `hg19` or `hg38` |
| `--groups` | Expression groups | `1,2,3,4` | No | `1,2,3,4,5` |
| `--output_prefix` | Output file prefix | None | Yes | `H1.3_profile` |
| `--output_dir` | Output directory | None | Yes | `./results` |
| `--window` | Window size (bp) | `3000` | No | `5000` |
| `--bin` | Bin size (bp) | `10` | No | `20` |

### 3. Stable Gene Analysis

```bash
# Basic usage
./run_stable_genes.sh

# Direct script usage
Rscript stable_genes.R \
  --parameters '{"file_path": "/path/to/data.rds", "n_groups": 10, "threshold": 0.15, "output_dir": "results"}'
```

#### Parameters

| Parameter | Description | Default | Required | Example |
|-----------|-------------|---------|----------|---------|
| `file_path` | Input RDS file | None | Yes | `All_Genes_TCGA_GTex_RSEM_TPM_converted.rds` |
| `n_groups` | Number of groups | `10` | No | `5` |
| `threshold` | Stability threshold | `0.15` | No | `0.2` |
| `output_dir` | Output directory | `results` | No | `./stable_genes` |

## Input/Output Formats

### Input Files

1. **RNA-seq Data**:
   - Format: Tab-separated values (TSV)
   - Required columns: Gene ID, Raw counts
   - Example: `GSE236538_RNASeq_H1Xsh_raw_read_counts_EnsemblID.txt`

2. **Histone Signal Data**:
   - Format: WIG or BigWig
   - Example: `GSM7578579_H1.3_WT.wig.gz`

3. **Gene Expression Data**:
   - Format: RDS (R data file)
   - Contains: TPM values for genes across samples
   - Example: `All_Genes_TCGA_GTex_RSEM_TPM_converted.rds`

### Output Files

1. **RNA-seq Quantile Results**:
   - Location: `quantile_results/`
   - Format: TSV files for each quantile group
   - Example: `RNASeq_quantile_1.tsv`

2. **TSS Profiles**:
   - Location: `results/`
   - Format: BigWig and matrix files
   - Example: `H1.3_profile_matrix.gz`

3. **Stable Gene Results**:
   - Location: `results/stable_genes/`
   - Format: RDS and CSV files
   - Example: `stable_genes_list.csv`

## SLURM Integration

The scripts are designed to work with SLURM job scheduler:

```bash
# Submit RNA-seq quantile analysis
sbatch run_rna_seq_quantile.sh

# Submit TSS profile generation (array job)
sbatch --array=0-9 run_deeptools_tss_profile_v2.sh

# Submit stable gene analysis
sbatch run_stable_genes.sh
```

### SLURM Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--time` | Job time limit | 4:00:00 |
| `--mem` | Memory per node | 20G |
| `--cpus-per-task` | CPU cores | 6 |
| `--partition` | Compute partition | normal |

## Recommended Workflow

1. **Prepare RNA-seq Data**:
   ```bash
   ./run_rna_seq_quantile.sh
   ```

2. **Generate TSS Profiles**:
   ```bash
   sbatch --array=0-9 run_deeptools_tss_profile_v2.sh
   ```

3. **Analyze Stable Genes**:
   ```bash
   ./run_stable_genes.sh
   ```

## Differences from Main Pipeline

1. **DeepTools Integration**:
   - Uses deepTools for profile generation
   - Better handling of large datasets
   - More efficient memory usage

2. **Additional Features**:
   - TPM normalization
   - Stable gene analysis
   - Parallel processing support

3. **Output Formats**:
   - Compatible with deepTools visualization
   - Additional statistical outputs
   - More detailed metadata

## Troubleshooting

1. **Memory Issues**:
   - Adjust `--mem` in SLURM parameters
   - Use `--bin` to increase bin size
   - Process fewer groups at a time

2. **Timeouts**:
   - Increase `--time` in SLURM parameters
   - Reduce window size with `--window`
   - Process smaller gene sets

3. **File Format Issues**:
   - Ensure correct file formats
   - Check file permissions
   - Verify genome versions match

## Contact

For questions or issues, please contact: chhetribsurya@gmail.com
