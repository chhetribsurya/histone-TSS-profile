#!/usr/bin/env Rscript

# Load required packages
library(optparse)
library(dplyr)
library(matrixStats)
library(jsonlite)

# Default input data directory and file
Data_dir <- '/data/baca/users/ss2976/GeneExp_Project/M2GEP_Shared_Repo/Data/Other_Data_Sources/'
DEFAULT_FILE_NAME <- 'All_Genes_TCGA_GTex_RSEM_TPM_converted.rds'
DEFAULT_FILE_PATH <- paste0(Data_dir, DEFAULT_FILE_NAME)

#' Load gene expression data from an RDS file
#'
#' @param file_path Character string specifying the path to the RDS file. 
#'                  Defaults to DEFAULT_FILE_PATH if not provided.
#' @return A data frame containing gene expression data with gene information 
#'         and sample expression values.
#' @details This function reads an RDS file containing gene expression data.
#'          The expected format includes columns for gene information 
#'          ('gene_id_no_version', 'gene_name', 'gene_id', 'chr', 'start', 
#'          'end', 'strand') and sample expression data (columns starting 
#'          with 'GTEX' for healthy samples or 'TCGA' for cancer samples).
#'          The function checks if the file exists and stops with an error 
#'          if it does not.
#' @examples
#' \dontrun{
#'   df <- load_data("/path/to/data.rds")
#' }
#' 
# Function to load gene expression data
load_data <- function(file_path = DEFAULT_FILE_PATH) {
    # Load gene expression data from an RDS file
    #
    # Parameters:
    #   file_path: str, path to the RDS file (default: DEFAULT_FILE_PATH)
    #
    # Returns:
    #   data.frame containing the gene expression data
    #
    # Expected Input Format:
    #   The input file should be an RDS file containing a data frame with:
    #   - Columns for gene information: 'gene_id_no_version', 'gene_name', 'gene_id', 'chr', 'start', 'end', 'strand'
    #   - Columns for sample expression data starting with 'GTEX' (healthy samples) or 'TCGA' (cancer samples)
    #   - The matrix should have genes as rows and samples as columns, with expression values (e.g., TPM) as the data
    if (!file.exists(file_path)) {
        stop(paste("Data file not found at:", file_path))
    }
    DF_TCGA_GTex <- readRDS(file_path)
    return(DF_TCGA_GTex)
}

#' Identify column types in the gene expression data frame
#'
#' @param df Data frame containing gene expression data.
#' @return A list with four elements: 
#'         - sample_cols: vector of all sample column names
#'         - healthy_cols: vector of healthy sample column names (GTEX)
#'         - cancer_cols: vector of cancer sample column names (TCGA)
#'         - info_cols: vector of gene information column names
#' @details This function categorizes columns in the data frame into sample 
#'          columns (starting with 'GTEX' or 'TCGA'), healthy columns 
#'          ('GTEX'), cancer columns ('TCGA'), and gene information columns.
#' @examples
#' \dontrun{
#'   cols <- identify_columns(df)
#' }
#' 
# Function to identify column types
identify_columns <- function(df) {
    # Identify sample, healthy, cancer, and information columns
    #
    # Parameters:
    #   df: data.frame, input data
    #
    # Returns:
    #   list with four elements: sample_cols, healthy_cols, cancer_cols, info_cols
    sample_cols <- grep("^(GTEX|TCGA)", names(df), value = TRUE)
    healthy_cols <- grep("^GTEX", names(df), value = TRUE)
    cancer_cols <- grep("^TCGA", names(df), value = TRUE)
    info_cols <- c("gene_id_no_version", "gene_name", "gene_id", "chr", "start", "end", "strand")
    list(sample_cols = sample_cols, 
         healthy_cols = healthy_cols, 
         cancer_cols = cancer_cols, 
         info_cols = info_cols)
}

#' Calculate statistical measures for gene expression data
#'
#' @param df Data frame containing gene expression data.
#' @param sample_cols Vector of column names for all samples.
#' @param healthy_cols Vector of column names for healthy samples.
#' @param cancer_cols Vector of column names for cancer samples.
#' @return Data frame with additional columns for averages, standard 
#'         deviations, and coefficient of variation (CV) across all samples, 
#'         healthy samples, and cancer samples.
#' @details This function computes row-wise means and standard deviations 
#'          for specified sample columns, and calculates the CV as SD/mean 
#'          for all samples combined.
#' @examples
#' \dontrun{
#'   df_stats <- calculate_statistics(df, cols$sample_cols, cols$healthy_cols, cols$cancer_cols)
#' }
#' 
# Function to calculate statistics
calculate_statistics <- function(df, sample_cols, healthy_cols, cancer_cols) {
    # Calculate averages, standard deviations, and CV for gene expression
    #
    # Parameters:
    #   df: data.frame, input data
    #   sample_cols: vector, columns for all samples
    #   healthy_cols: vector, columns for healthy samples
    #   cancer_cols: vector, columns for cancer samples
    #
    # Returns:
    #   data.frame with additional columns for statistics
    df$AVERAGE_All <- rowMeans(df[, sample_cols], na.rm = TRUE)
    df$SD_All <- rowSds(as.matrix(df[, sample_cols]), na.rm = TRUE)
    df$AVERAGE_GTEX <- rowMeans(df[, healthy_cols], na.rm = TRUE)
    df$SD_GTEX <- rowSds(as.matrix(df[, healthy_cols]), na.rm = TRUE)
    df$AVERAGE_TCGA <- rowMeans(df[, cancer_cols], na.rm = TRUE)
    df$SD_TCGA <- rowSds(as.matrix(df[, cancer_cols]), na.rm = TRUE)
    df$CV_All <- df$SD_All / abs(df$AVERAGE_All)
    df
}

#' Group genes into bins based on average expression levels
#'
#' @param df Data frame with an AVERAGE_All column.
#' @param n_groups Integer specifying the number of groups (default: 10).
#' @return Data frame with an additional Group_Num column indicating the 
#'         expression level group for each gene.
#' @details This function sorts genes by their average expression across all 
#'          samples and assigns them to equal-sized bins using the ntile 
#'          function from dplyr.
#' @examples
#' \dontrun{
#'   df_grouped <- group_genes(df, n_groups = 10)
#' }
#' 
# Function to group genes by expression level
group_genes <- function(df, n_groups = 10) {
    # Group genes into equal-sized bins based on average expression
    #
    # Parameters:
    #   df: data.frame, input data with AVERAGE_All column
    #   n_groups: int, number of groups (default: 10)
    #
    # Returns:
    #   data.frame with an additional Group_Num column
    df <- df %>% arrange(AVERAGE_All)
    df$Group_Num <- ntile(df$AVERAGE_All, n_groups)
    df
}

#' Select stable genes based on coefficient of variation
#'
#' @param df Data frame with Group_Num and CV_All columns.
#' @param threshold Float between 0 and 1 specifying the proportion of genes 
#'                  to select per group (default: 0.15).
#' @return Data frame containing the selected stable genes with the lowest 
#'         CV within each group.
#' @details This function groups genes by their assigned Group_Num, sorts 
#'          them by CV within each group, and selects the top proportion 
#'          (threshold) of genes with the lowest CV.
#' @examples
#' \dontrun{
#'   stable_df <- select_stable_genes(df, threshold = 0.15)
#' }
#' 
# Function to select stable genes
select_stable_genes <- function(df, threshold = 0.05) {
    # Select genes with the lowest CV within each group
    #
    # Parameters:
    #   df: data.frame, input data with Group_Num and CV_All columns
    #   threshold: float, proportion of genes to select per group (default: 0.15)
    #
    # Returns:
    #   data.frame containing the selected stable genes
    stable_genes_df <- df %>%
        group_by(Group_Num) %>%
        arrange(CV_All) %>%
        slice_head(prop = threshold) %>%
        ungroup()
    stable_genes_df
}

#' Main function to orchestrate stable gene selection
#'
#' @param parameters List containing:
#'                   - file_path: path to input RDS file (optional)
#'                   - n_groups: number of expression groups
#'                   - threshold: proportion of stable genes to select
#'                   - output_dir: directory to save results
#' @return None. Outputs results to console and saves stable genes to a CSV 
#'         file in the specified output directory.
#' @details This function coordinates the entire analysis pipeline: loading 
#'          data, calculating statistics, grouping genes, selecting stable 
#'          genes, and saving results. It includes error handling for file 
#'          existence and directory creation, and prints an interpretability 
#'          statement.
#' @examples
#' \dontrun{
#'   parameters <- list(file_path = "data.rds", n_groups = 10, threshold = 0.50, output_dir = "results")
#'   main(parameters)
#' }
#' 
# Main function to orchestrate the analysis
main <- function(parameters) {
    # Orchestrate the process of selecting stable genes
    #
    # Parameters:
    #   parameters: list with file_path, n_groups, threshold, output_dir
    #
    # The output_dir parameter specifies where the stable genes will be saved
    file_path <- if (!is.null(parameters$file_path)) parameters$file_path else DEFAULT_FILE_PATH
    n_groups <- parameters$n_groups
    threshold <- parameters$threshold
    output_dir <- parameters$output_dir
    
    # Print parameters used for processing
    cat("Parameters used:\n")
    cat(sprintf("  file_path: %s\n", file_path))
    cat(sprintf("  n_groups: %d\n", n_groups))
    cat(sprintf("  threshold: %.2f\n", threshold))
    cat(sprintf("  output_dir: %s\n", output_dir))
    
    # Load data and process it
    df <- load_data(file_path)
    cols <- identify_columns(df)
    df <- calculate_statistics(df, cols$sample_cols, cols$healthy_cols, cols$cancer_cols)
    df <- group_genes(df, n_groups)
    stable_genes_df <- select_stable_genes(df, threshold)
    select_cols_stable_genes_df <- stable_genes_df[, c("gene_id_no_version", "gene_name", "AVERAGE_All", "SD_All", "CV_All", "Group_Num")]
    
    # Output results
    stable_genes <- stable_genes_df$gene_id_no_version
    cat(sprintf("Selected %d stable genes.\n", length(stable_genes)))
    
    # Print a preview of the final result dataframe (full columns)
    cat("\nFirst few rows of the stable genes dataframe (full columns):\n")
    print(head(stable_genes_df))
    
    # Print a preview of the final result dataframe
    cat("\nFirst few rows of the stable genes dataframe (selected columns):\n")
    print(head(stable_genes_df[, c("gene_id_no_version", "gene_name", "AVERAGE_All", "SD_All", "CV_All", "Group_Num")]))
    
    # Print dimensions and column names
    cat(sprintf("\nDimensions of stable genes dataframe: %d rows, %d columns\n", nrow(stable_genes_df), ncol(stable_genes_df)))
    cat("Columns in the stable genes dataframe:\n")
    print(names(stable_genes_df))
    
    # Print number of stable genes per group
    cat("\nNumber of stable genes per group:\n")
    print(table(stable_genes_df$Group_Num))
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Save results to file in the specified output directory
    output_file1 <- file.path(output_dir, paste0("stable_genes_", threshold, ".tsv"))
    output_file2 <- file.path(output_dir, paste0("stable_genes_select_cols_", threshold, ".tsv"))
    write.table(stable_genes_df, output_file1, row.names = FALSE, quote = FALSE, sep="\t")
    write.table(select_cols_stable_genes_df, output_file2, row.names = FALSE, quote = FALSE, sep="\t")
    cat(sprintf("Results saved to %s\n", output_file2))
    
    # Explain interpretability
    cat("\nInterpretability of Results:\n")
    cat("The selected stable genes are those with the lowest coefficient of variation (CV) within each expression level group.\n")
    cat("These genes exhibit consistent expression levels relative to their mean expression, making them reliable for normalization or as reference genes.\n")
}

# Define command-line options
option_list <- list(
    make_option(c("--parameters"), type = "character", 
                default = paste0('{"file_path": "', DEFAULT_FILE_PATH, '", "n_groups": 10, "threshold": 0.05, "output_dir": "results"}'), 
                help = "JSON string with parameters: file_path (str), n_groups (int), threshold (float), output_dir (str)")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Parse JSON parameters with error handling
tryCatch({
    parameters <- fromJSON(opt$parameters)
}, error = function(e) {
    stop("Invalid JSON in --parameters: ", e$message)
})

# Execute the main function
main(parameters)
