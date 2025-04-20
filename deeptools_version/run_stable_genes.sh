#!/usr/bin/bash

# Default input data directory and file
DATA_DIR='/data/baca/users/ss2976/GeneExp_Project/M2GEP_Shared_Repo/Data/Other_Data_Sources/'
DEFAULT_FILE_NAME='All_Genes_TCGA_GTex_RSEM_TPM_converted.rds'
DEFAULT_FILE_PATH="${DATA_DIR}/${DEFAULT_FILE_NAME}"

#Rscript stable_genes.R --parameters '{"file_path": "/path/to/your/data.rds", "n_groups": 10, "threshold": 0.15, "output_dir": "results"}'
Rscript stable_genes.R --parameters '{"file_path": $DEFAULT_FILE_PATH, "n_groups": 10, "threshold": 0.15, "output_dir": "results"}'
