# Histone TSS Profile Analysis Pipeline

## Overview

This repository contains a comprehensive pipeline for analyzing histone variant distribution around transcription start sites (TSS) and enhancer RNAs (eRNAs). The pipeline processes genomic data to generate TSS profiles, analyze histone distribution patterns, and create smoothed visualizations.

### Biological Context

Histone variants play crucial roles in gene regulation and chromatin organization. Different histone variants, such as H1.0, H1.2, H1.3, H1.4, H1.5, and H1.X, have distinct distribution patterns around transcription start sites (TSS) and enhancer regions. Understanding these patterns provides insights into:
- Gene regulation mechanisms
- Chromatin structure and dynamics
- Transcription factor binding
- Enhancer-promoter interactions

This pipeline enables systematic analysis of histone variant distribution patterns across different gene expression groups and enhancer regions, providing valuable insights into chromatin organization and gene regulation.

## Pipeline Overview

The analysis pipeline consists of three main steps:

1. **TSS BED Generation** (`step1_*` scripts)
   - Generates TSS BED files from RefSeq or GENCODE annotations
   - Groups genes by expression levels
   - Creates input files for downstream analysis

2. **Histone Distribution Analysis** (`step2_*` scripts)
   - Analyzes histone variant distribution around TSS
   - Processes WIG/BigWig files
   - Generates profile plots and statistical analyses

3. **Profile Smoothing** (`step3_*` scripts)
   - Applies multiple smoothing methods to TSS profiles
   - Compares different smoothing parameters
   - Generates publication-quality visualizations

## Input/Output Formats

### Input Files

1. **Gene Expression Data**:
   - Format: Tab-separated values (TSV)
   - Required columns: Gene ID, Expression Value
   - Location: `./data/T47D_ExpressionGroups_RNASeq/`

2. **Histone Signal Data**:
   - Format: WIG or BigWig
   - Example: `GSM5076926_H1.0_WT.wig`, `GSM5076927_H1.2_WT.wig`
   - Location: `./data/`

3. **Genome Annotation**:
   - Format: BED
   - Source: RefSeq or GENCODE annotations
   - Genome versions supported: hg19, hg38

### Output Files

1. **TSS BED Files**:
   - Format: BED6+ (6 standard BED columns + additional metadata)
   - Columns: chrom, start, end, name, score, strand, group
   - Location: `refseq_tss_results/` or `gencode_tss_bed/`

2. **Profile Data**:
   - Format: CSV
   - Contains: Position, Signal Value, Group
   - Location: `resultsFINAL/results_H1.*/h1_tss_profile_profiles.csv`

3. **Visualization Files**:
   - Format: PNG, PDF
   - Types: Profile plots, smoothed profiles
   - Location: `resultsFINAL/results_H1.*/` and `smoothing_comparison/`

## Directory Structure

```
histone-TSS-profile/
├── data/                    # Input data directory
│   ├── T47D_ExpressionGroups_RNASeq/  # Gene expression data
│   └── histone_signal/      # Histone signal files (WIG/BigWig)
├── deeptools_version/       # Alternative implementation using deepTools
├── stable_genes_version/    # Version using stable gene sets
├── scripts/                 # Main analysis scripts
├── resultsFINAL/           # Output directory for final results
└── LICENSE                 # License information
```

## Command-Line Usage

### Step 1: TSS BED Generation

```bash
# For RefSeq annotations
./step1_run_refseq-tss-bed-generator.sh \
  --expr_dir ./data/T47D_ExpressionGroups_RNASeq \
  --output_dir refseq_tss_results \
  --genome hg19

# For eRNA analysis
./step1_run_erna-tss-bed-generator.sh \
  --input ./data/enhancer_rna_files/ \
  --output_dir eRNA_tss_results
```

### Step 2: Histone Distribution Analysis

```bash
./step2_run_h1_tss_profile_v2.sh \
  --tss ./refseq_tss_results/refseq_tss_with_groups.bed \
  --input ./data/GSM5076926_H1.0_WT.wig \
  --window 3000 \
  --bin 10 \
  --processes 6 \
  --groups "1,2,3,4"
```

### Step 3: Profile Smoothing

```bash
./step3_run_smooth_tss_profiles.sh \
  --input ./resultsFINAL/results_H1.0/h1_tss_profile_profiles.csv \
  --output-dir ./resultsFINAL/results_H1.0/smoothing_comparison \
  --method savgol \
  --window 21 \
  --poly 3
```

## Main Scripts

### Step 1: TSS BED Generation
- `step1_run_refseq-tss-bed-generator.sh`: Generates TSS BED files from RefSeq annotations
- `step1_run_erna-tss-bed-generator.sh`: Generates TSS BED files for eRNAs
- `refseq-tss-bed-generator.py`: Python script for RefSeq TSS processing
- `gencode_tss_bed_generator_v2.py`: Python script for GENCODE TSS processing

### Step 2: Histone Distribution Analysis
- `step2_run_h1_tss_profile_v2.sh`: Main script for histone distribution analysis
- `h1_tss_profile.py`: Python script for processing histone data

### Step 3: Profile Smoothing
- `step3_run_smooth_tss_profiles.sh`: Script for profile smoothing
- `smooth_tss_profiles.py`: Python script implementing smoothing methods

## Usage

### Prerequisites
- Python 3.x
- Required Python packages (install via conda):
  ```bash
  conda create -n genometools_env python=3.x
  conda activate genometools_env
  conda install -c bioconda bedtools pybedtools numpy pandas matplotlib scipy
  ```

### Running the Pipeline

1. **Generate TSS BED files**:
   ```bash
   ./step1_run_refseq-tss-bed-generator.sh
   # or
   ./step1_run_erna-tss-bed-generator.sh
   ```

2. **Analyze histone distribution**:
   ```bash
   ./step2_run_h1_tss_profile_v2.sh
   ```

3. **Smooth profiles**:
   ```bash
   ./step3_run_smooth_tss_profiles.sh
   ```

## Configuration Parameters

### Step 1: TSS BED Generation Parameters

| Parameter | Description | Default | Required | Example |
|-----------|-------------|---------|----------|---------|
| `--expr_dir` | Directory containing expression data | None | Yes | `./data/T47D_ExpressionGroups_RNASeq` |
| `--output_dir` | Output directory for BED files | `refseq_tss_results` | No | `./results/tss_bed` |
| `--genome` | Genome version | `hg19` | No | `hg19` or `hg38` |
| `--groups` | Expression groups to analyze | `1,2,3,4` | No | `1,2,3,4,5` |
| `--min_expr` | Minimum expression value | `0` | No | `1.0` |
| `--max_expr` | Maximum expression value | `None` | No | `100.0` |
| `--group_size` | Number of genes per group | `None` | No | `1000` |

### Step 2: Histone Distribution Analysis Parameters

| Parameter | Description | Default | Required | Example |
|-----------|-------------|---------|----------|---------|
| `--tss` | TSS BED file | None | Yes | `./refseq_tss_results/refseq_tss_with_groups.bed` |
| `--input` | Input WIG/BigWig file | None | Yes | `./data/GSM5076926_H1.0_WT.wig` |
| `--window` | Window size around TSS (bp) | `3000` | No | `5000` |
| `--bin` | Bin size for signal averaging (bp) | `10` | No | `20` |
| `--processes` | Number of CPU cores | `1` | No | `6` |
| `--groups` | Expression groups to analyze | `1,2,3,4` | No | `1,2,3,4,5` |
| `--convert` | Convert WIG to BigWig | `False` | No | `True` |
| `--chrom_sizes` | Chromosome sizes file | None | If `--convert` | `hg19.chrom.sizes` |
| `--sample` | Number of TSS to sample | `None` | No | `10000` |
| `--format` | Input file format | `wig` | No | `bigwig` |
| `--abs` | Use absolute values | `False` | No | `True` |
| `--log2` | Log2 transform signal | `False` | No | `True` |

### Step 3: Profile Smoothing Parameters

| Parameter | Description | Default | Required | Example |
|-----------|-------------|---------|----------|---------|
| `--input` | Input profile CSV file | None | Yes | `./resultsFINAL/results_H1.0/h1_tss_profile_profiles.csv` |
| `--output-dir` | Output directory | None | Yes | `./resultsFINAL/results_H1.0/smoothing_comparison` |
| `--method` | Smoothing method | `savgol` | No | `gaussian` or `moving_avg` |
| `--window` | Window size for smoothing | `21` | No | `51` |
| `--sigma` | Sigma for Gaussian smoothing | `1` | If `method=gaussian` | `2` |
| `--poly` | Polynomial order for Savitzky-Golay | `3` | If `method=savgol` | `4` |
| `--prefix` | Output file prefix | `h1_smoothed` | No | `h1_gaussian_smoothed` |

### Common Parameters Across Steps

| Parameter | Description | Default | Required | Example |
|-----------|-------------|---------|----------|---------|
| `--help` | Show help message | N/A | No | N/A |
| `--version` | Show version information | N/A | No | N/A |
| `--verbose` | Enable verbose output | `False` | No | `True` |
| `--log` | Log file path | `None` | No | `./analysis.log` |

### Expression Group Definitions

| Group | Description | Expression Range | Typical Usage |
|-------|-------------|------------------|---------------|
| 0 | No expression | < 1.0 | Usually excluded |
| 1 | Low expression | 1.0 - 10.0 | Baseline group |
| 2 | Medium expression | 10.0 - 100.0 | Intermediate group |
| 3 | High expression | 100.0 - 1000.0 | Highly expressed |
| 4 | Very high expression | > 1000.0 | Maximum expression |

### Smoothing Method Comparison

| Method | Description | Best Use Case | Parameters |
|--------|-------------|---------------|------------|
| `savgol` | Savitzky-Golay filter | General purpose smoothing | window, poly |
| `gaussian` | Gaussian smoothing | Noise reduction | window, sigma |
| `moving_avg` | Moving average | Simple smoothing | window |

### Recommended Parameter Combinations

1. **Basic Analysis**:
   ```bash
   --window 3000 --bin 10 --processes 4 --groups 1,2,3,4
   ```

2. **High-Resolution Analysis**:
   ```bash
   --window 5000 --bin 5 --processes 8 --groups 1,2,3,4,5
   ```

3. **Quick Analysis**:
   ```bash
   --window 2000 --bin 20 --processes 2 --groups 1,2,3
   ```

4. **Smoothing Parameters**:
   - Savitzky-Golay: `--window 21 --poly 3`
   - Gaussian: `--window 31 --sigma 2`
   - Moving Average: `--window 51`

## Output Files

- **TSS BED Generation**:
  - `refseq_tss_results/`: Contains generated TSS BED files
  - `gencode_tss_bed/`: Contains GENCODE-based TSS BED files

- **Histone Distribution**:
  - `resultsFINAL/results_H1.*/`: Contains profile plots and analysis results
  - `h1_tss_profile_profiles.csv`: Raw profile data

- **Profile Smoothing**:
  - `smoothing_comparison/`: Contains smoothed profiles using different methods
  - Subdirectories for each smoothing method (savgol, gaussian, moving_avg)

## Biological Interpretation

The pipeline generates several types of analyses that provide biological insights:

1. **TSS Distribution Patterns**:
   - Shows how different histone variants are distributed around TSS
   - Reveals variant-specific binding preferences
   - Identifies potential regulatory roles

2. **Expression Group Analysis**:
   - Compares histone variant distribution across expression levels
   - Identifies variant-specific patterns in gene regulation
   - Reveals potential roles in transcriptional control

3. **Enhancer Analysis**:
   - Analyzes histone variant distribution around eRNAs
   - Provides insights into enhancer regulation
   - Reveals variant-specific roles in enhancer function

## License

This project is licensed under the terms specified in the LICENSE file.

## Contact

For questions, improvements, or issues, please contact: chhetribsurya@gmail.com
