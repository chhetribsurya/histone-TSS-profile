#!/usr/bin/env Rscript

# Load required packages
library(optparse)
library(dplyr)
library(matrixStats)
library(jsonlite)

#' Load gene expression data from an RDS file
#'
#' @param file_path Character string specifying the path to the RDS file.
#' @return A data frame containing gene expression data with gene information 
#'         and sample expression values.
#' 
# Function to load gene expression data
load_data <- function(file_path) {
    # Check if file exists
    if (!file.exists(file_path)) {
        stop(paste("Data file not found at:", file_path))
    }
    
    # Load the RDS file
    DF <- readRDS(file_path)
    return(DF)
}

#' Identify column types in the gene expression data frame
#'
#' @param df Data frame containing gene expression data.
#' @return A list with four elements: 
#'         - sample_cols: vector of all sample column names
#'         - healthy_cols: vector of healthy sample column names (GTEX)
#'         - cancer_cols: vector of cancer sample column names (TCGA)
#'         - info_cols: vector of gene information column names
#' 
# Function to identify column types
identify_columns <- function(df) {
    # Identify sample, healthy, cancer, and information columns
    sample_cols <- grep("^(GTEX|TCGA)", names(df), value = TRUE)
    healthy_cols <- grep("^GTEX", names(df), value = TRUE)
    cancer_cols <- grep("^TCGA", names(df), value = TRUE)
    info_cols <- c("gene_id_no_version", "gene_name", "gene_id", "chr", "start", "end", "strand")
    list(sample_cols = sample_cols, 
         healthy_cols = healthy_cols, 
         cancer_cols = cancer_cols, 
         info_cols = info_cols)
}

#' Calculate average expression for gene expression data
#'
#' @param df Data frame containing gene expression data.
#' @param sample_cols Vector of column names for all samples.
#' @return Data frame with additional columns for averages across all samples.
#' 
# Function to calculate average expression
calculate_average_expression <- function(df, sample_cols) {
    # Calculate average expression for all samples
    df$AVERAGE_All <- rowMeans(df[, sample_cols], na.rm = TRUE)
    df
}

#' Group genes into quantiles based on average expression levels
#'
#' @param df Data frame with an AVERAGE_All column.
#' @param n_groups Integer specifying the number of quantile groups.
#' @return Data frame with an additional Group_Num column indicating the 
#'         expression level group for each gene.
#' 
# Function to group genes by expression level into quantiles
group_genes_by_quantile <- function(df, n_groups) {
    # Group genes into equal-sized bins based on average expression
    df <- df %>% arrange(AVERAGE_All)
    df$Group_Num <- ntile(df$AVERAGE_All, n_groups)
    df
}

#' Save gene IDs for each quantile group to separate files
#'
#' @param df Data frame with Group_Num and gene_id_no_version columns.
#' @param output_dir Directory where output files will be saved.
#' @param prefix Prefix for the output files.
#' @return None. Files are saved to the specified directory.
#' 
# Function to save gene IDs for each quantile group
save_quantile_groups <- function(df, output_dir, prefix, id_column = "gene_id_no_version") {
    # Make sure the output directory exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Get all unique group numbers
    groups <- unique(df$Group_Num)
    
    # Sort groups numerically
    groups <- sort(groups)
    
    # For each group, save a file with the gene IDs
    for (group in groups) {
        # Get gene IDs for this group
        group_genes <- df %>% 
            filter(Group_Num == group) %>% 
            pull(!!sym(id_column))
        
        # Create the output file path
        output_file <- file.path(output_dir, paste0(prefix, "_ENS_", group - 1, ".txt"))
        
        # Save the gene IDs to a file (one per line, no header)
        writeLines(group_genes, output_file)
        
        cat(sprintf("Saved %d genes for group %d to %s\n", length(group_genes), group, output_file))
    }
    
    # Save all genes in one file with their group numbers
    all_genes_file <- file.path(output_dir, paste0(prefix, "_all_groups.tsv"))
    write.table(df[, c(id_column, "Group_Num", "AVERAGE_All")], 
                all_genes_file, 
                row.names = FALSE, 
                quote = FALSE,
                sep = "\t")
    
    cat(sprintf("Saved all genes with group numbers to %s\n", all_genes_file))
}

#' Main function to orchestrate gene quantile grouping
#'
#' @param file_path Path to input RDS file.
#' @param output_dir Directory to save results.
#' @param n_groups Number of quantile groups.
#' @param prefix Prefix for output filenames.
#' @param id_column Column name for gene identifiers.
#' @return None. Outputs results to console and saves gene groups to files.
#' 
# Main function to orchestrate the analysis
main <- function(file_path, output_dir, n_groups, prefix, id_column) {
    # Print parameters used for processing
    cat("Parameters used:\n")
    cat(sprintf("  file_path: %s\n", file_path))
    cat(sprintf("  output_dir: %s\n", output_dir))
    cat(sprintf("  n_groups: %d\n", n_groups))
    cat(sprintf("  prefix: %s\n", prefix))
    cat(sprintf("  id_column: %s\n", id_column))
    
    # Load data
    cat("\nLoading data...\n")
    df <- load_data(file_path)
    cat(sprintf("Loaded data with %d rows and %d columns\n", nrow(df), ncol(df)))
    
    # Identify column types
    cat("\nIdentifying column types...\n")
    cols <- identify_columns(df)
    cat(sprintf("Found %d sample columns\n", length(cols$sample_cols)))
    
    # Calculate average expression
    cat("\nCalculating average expression...\n")
    df <- calculate_average_expression(df, cols$sample_cols)
    
    # Group genes by expression quantiles
    cat("\nGrouping genes by expression quantiles...\n")
    df <- group_genes_by_quantile(df, n_groups)
    
    # Report number of genes per group
    cat("\nNumber of genes per expression quantile:\n")
    group_counts <- table(df$Group_Num)
    print(group_counts)
    
    # Save genes for each quantile group
    cat("\nSaving gene lists for each quantile group...\n")
    save_quantile_groups(df, output_dir, prefix, id_column)
    
    cat("\nProcessing complete!\n")
}

# Define command-line options
option_list <- list(
    make_option(c("--input"), type = "character", 
                help = "Path to input RDS file"),
    make_option(c("--output-dir"), type = "character", default = "results",
                help = "Directory to save output files [default: %default]"),
    make_option(c("--num-groups"), type = "integer", default = 10,
                help = "Number of quantile groups [default: %default]"),
    make_option(c("--prefix"), type = "character", default = "quantile",
                help = "Prefix for output files [default: %default]"),
    make_option(c("--id-column"), type = "character", default = "gene_id_no_version",
                help = "Column name for gene IDs [default: %default]")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$input)) {
    stop("Input file path is required. Use --input option.")
}

# Execute the main function
main(
    file_path = opt$input,
    output_dir = opt$`output-dir`,
    n_groups = opt$`num-groups`,
    prefix = opt$prefix,
    id_column = opt$`id-column`
)
