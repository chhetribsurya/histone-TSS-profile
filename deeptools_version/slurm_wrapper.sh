#!/usr/bin/bash

# This script counts the files and submits the SLURM array job with the correct range

# Directory with the gene files
GENE_FILES_DIR="/data/baca/users/sc1238/datasets/scripts/for_alexis/quantile_results_new_normalized/equal_bins"

# Count the actual number of files
FILE_COUNT=$(ls ${GENE_FILES_DIR}/RNASeq*.tsv | wc -l)

if [ $FILE_COUNT -eq 0 ]; then
    echo "Error: No matching files found!"
    exit 1
fi

# Calculate the maximum array index (0-based indexing)
MAX_INDEX=$((FILE_COUNT - 1))

echo "Found $FILE_COUNT files. Will run array jobs with indices 0-$MAX_INDEX"

# Submit the job with the dynamically determined array range
#sbatch --array=0-${MAX_INDEX}%2 run_deeptools_tss_profile_v2.sh
sbatch --array=0-${MAX_INDEX} run_deeptools_tss_profile_v2.sh
