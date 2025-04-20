#!/usr/bin/env Rscript

# Load required packages
library(optparse)
library(dplyr)
library(matrixStats)

#' Load gene expression data from an RDS file
#'
#' @param file_path Character string specifying the path to the RDS file.
#' @return A data frame containing gene expression data
load_data <- function(file_path) {
    if (!file.exists(file_path)) {
        stop(paste("Data file not found at:", file_path))
    }
    return(readRDS(file_path))
}

#' Identify column types in the gene expression data frame
#'
#' @param df Data frame containing gene expression data.
#' @return A list with column type categories
identify_columns <- function(df) {
    sample_cols <- grep("^(GTEX|TCGA)", names(df), value = TRUE)
    healthy_cols <- grep("^GTEX", names(df), value = TRUE)
    cancer_cols <- grep("^TCGA", names(df), value = TRUE)
    info_cols <- c("gene_id_no_version", "gene_name", "gene_id", "chr", "start", "end", "strand")
    
    return(list(
        sample_cols = sample_cols, 
        healthy_cols = healthy_cols, 
        cancer_cols = cancer_cols, 
        info_cols = info_cols
    ))
}

#' Calculate average expression for genes
#'
#' @param df Data frame containing gene expression data
#' @param sample_cols Vector of column names for samples
#' @return Data frame with an AVERAGE_All column added
calculate_average_expression <- function(df, sample_cols) {
    df$AVERAGE_All <- rowMeans(df[, sample_cols], na.rm = TRUE)
    return(df)
}

#' Group genes into quantiles based on average expression
#'
#' @param df Data frame with an AVERAGE_All column
#' @param n_groups Integer specifying the number of quantile groups
#' @param include_non_expressed Logical, create a separate group for non-expressed genes
#' @param non_expressed_threshold Numeric, threshold for non-expressed genes
#' @return Data frame with Group_Num column added
group_genes_by_quantile <- function(df, n_groups, include_non_expressed = FALSE, non_expressed_threshold = 0.01) {
    if (include_non_expressed) {
        # Identify non-expressed genes
        non_expressed <- df$AVERAGE_All < non_expressed_threshold
        
        # Initialize all genes as group 0
        df$Group_Num <- 0
        
        # Only process expressed genes if there are any
        if (sum(!non_expressed) > 0) {
            # Calculate quantiles only for expressed genes
            expressed_df <- df[!non_expressed, ]
            expressed_df$Group_Num <- ntile(expressed_df$AVERAGE_All, n_groups)
            
            # Update the groups in the original dataframe
            df$Group_Num[!non_expressed] <- expressed_df$Group_Num
            
            # Log counts
            cat(sprintf("Non-expressed genes (Group 0): %d\n", sum(non_expressed)))
            cat(sprintf("Expressed genes (Groups 1-%d): %d\n", n_groups, sum(!non_expressed)))
        } else {
            cat("Warning: All genes are below the non-expressed threshold\n")
        }
    } else {
        # Standard quantile grouping without separate non-expressed group
        df$Group_Num <- ntile(df$AVERAGE_All, n_groups)
    }
    
    return(df)
}

#' Save gene IDs for each quantile group
#'
#' @param df Data frame with Group_Num and ID columns
#' @param output_dir Directory for output files
#' @param prefix Prefix for output filenames
#' @param id_column Column name for gene IDs
#' @return None (files are saved)
save_quantile_groups <- function(df, output_dir, prefix, id_column = "gene_id_no_version") {
    # Create output directory if needed
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Get all unique group numbers
    groups <- sort(unique(df$Group_Num))
    
    # For each group, save a file with gene IDs
    for (group in groups) {
        # Get gene IDs for this group
        group_genes <- df %>% 
            filter(Group_Num == group) %>% 
            pull(!!sym(id_column))
        
        # Use the group number directly in the filename
        file_index <- group
        
        # Create output file path
        output_file <- file.path(output_dir, paste0(prefix, "_ENS_", file_index, ".txt"))
        
        # Save gene IDs to file (one per line, no header)
        writeLines(group_genes, output_file)
        
        cat(sprintf("Saved %d genes for group %d to %s\n", 
                    length(group_genes), group, output_file))
    }
    
    # Save combined file with all genes and groups
    all_genes_file <- file.path(output_dir, paste0(prefix, "_all_groups.tsv"))
    write.table(df[, c(id_column, "Group_Num", "AVERAGE_All")], 
                all_genes_file, 
                row.names = FALSE, 
                quote = FALSE,
                sep = "\t")
    
    cat(sprintf("Saved all genes with group numbers to %s\n", all_genes_file))
}

#' Main function to orchestrate the analysis
#'
#' @param file_path Path to input RDS file
#' @param output_dir Directory to save results
#' @param n_groups Number of quantile groups
#' @param prefix Prefix for output filenames
#' @param id_column Column name for gene IDs
#' @param include_non_expressed Include separate group for non-expressed genes
#' @param non_expressed_threshold Expression threshold for non-expressed genes
#' @return None (processing is done and files saved)
main <- function(file_path, output_dir, n_groups, prefix, id_column, 
                include_non_expressed = FALSE, non_expressed_threshold = 0.01) {
    # Print parameters
    cat("Parameters used:\n")
    cat(sprintf("  file_path: %s\n", file_path))
    cat(sprintf("  output_dir: %s\n", output_dir))
    cat(sprintf("  n_groups: %d\n", n_groups))
    cat(sprintf("  prefix: %s\n", prefix))
    cat(sprintf("  id_column: %s\n", id_column))
    cat(sprintf("  include_non_expressed: %s\n", ifelse(include_non_expressed, "TRUE", "FALSE")))
    if (include_non_expressed) {
        cat(sprintf("  non_expressed_threshold: %f\n", non_expressed_threshold))
    }
    
    # Run the analysis
    cat("\nLoading data...\n")
    df <- load_data(file_path)
    cat(sprintf("Loaded data with %d rows and %d columns\n", nrow(df), ncol(df)))
    
    cat("\nIdentifying column types...\n")
    cols <- identify_columns(df)
    cat(sprintf("Found %d sample columns\n", length(cols$sample_cols)))
    
    cat("\nCalculating average expression...\n")
    df <- calculate_average_expression(df, cols$sample_cols)
    
    cat("\nGrouping genes by expression quantiles...\n")
    df <- group_genes_by_quantile(df, n_groups, include_non_expressed, non_expressed_threshold)
    
    cat("\nNumber of genes per expression quantile:\n")
    print(table(df$Group_Num))
    
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
                help = "Column name for gene IDs [default: %default]"),
    make_option(c("--include-non-expressed"), action = "store_true", default = FALSE,
                help = "Include a separate group (0) for non-expressed genes [default: %default]"),
    make_option(c("--non-expressed-threshold"), type = "double", default = 0.01,
                help = "Expression threshold below which genes are considered non-expressed [default: %default]")
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
    id_column = opt$`id-column`,
    include_non_expressed = opt$`include-non-expressed`,
    non_expressed_threshold = opt$`non-expressed-threshold`
)
