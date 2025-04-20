#!/usr/bin/env python3
"""
RNA-seq Quantile Analysis Script with Integrated GENCODE Gene Length Extraction

This script processes RNA-seq data to group genes into 10 quantiles based on expression 
levels for each condition. It can normalize raw counts to TPM and automatically 
download/extract gene lengths from GENCODE v37 if needed.
"""

import pandas as pd
import numpy as np
import argparse
import os
import sys
import re
import gzip
import subprocess
import time
from urllib.request import urlretrieve

# GENCODE v37 URL
GENCODE_V37_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz"


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process RNA-seq data and group genes into quantiles.')
    parser.add_argument('--input', '-i', required=True, help='Input RNA-seq file path')
    parser.add_argument('--output_dir', '-o', required=True, help='Output directory for results')
    parser.add_argument('--quantiles', '-q', type=int, default=10, help='Number of quantiles (default: 10)')
    parser.add_argument('--id_col', default='ID', help='Column name for gene IDs (default: ID)')
    parser.add_argument('--min_expr', type=float, default=1.0, 
                        help='Minimum expression threshold to consider a gene expressed (default: 1.0)')
    parser.add_argument('--gene_name_file', '-g', 
                        help='Optional file mapping gene IDs to gene names (tab-separated)')
    parser.add_argument('--equal_bins', action='store_true', 
                        help='Also create strictly equal-sized bins for expressed genes (saved in subdirectory)')
    parser.add_argument('--normalize', action='store_true',
                        help='Normalize raw counts to TPM values before quantile assignment')
    parser.add_argument('--normalize_matrix_only', action='store_true',
                        help='Only perform TPM normalization of the full matrix and skip quantile analysis')
    parser.add_argument('--gene_length_file', 
                        help='Gene length file for TPM normalization (if not provided, will extract from GENCODE v37)')
    parser.add_argument('--gtf_file',
                        help='Optional GENCODE GTF file for gene length extraction (will download if not provided)')
    parser.add_argument('--no_download', action='store_true',
                        help='Do not download GENCODE GTF file automatically (will error if needed for normalization)')
    return parser.parse_args()


def load_rnaseq_data(file_path, id_col='ID'):
    """
    Load RNA-seq data from file.
    
    Args:
        file_path: Path to RNA-seq data file
        id_col: Column name containing gene IDs
        
    Returns:
        Pandas DataFrame with RNA-seq data
    """
    print(f"Loading RNA-seq data from {file_path}...")
    df = pd.read_csv(file_path, sep='\t', index_col=None)
    
    # Clean up gene IDs (remove quotes if present)
    if df[id_col].dtype == 'object':
        df[id_col] = df[id_col].str.replace('"', '')
    
    return df


def load_gene_names(file_path):
    """
    Load gene ID to gene name mapping.
    
    Args:
        file_path: Path to gene mapping file
        
    Returns:
        Dictionary mapping gene IDs to gene names
    """
    print(f"Loading gene name mappings from {file_path}...")
    mapping_df = pd.read_csv(file_path, sep='\t')
    
    # Assuming the file has columns 'gene_id' and 'gene_name'
    # Adjust column names as needed based on your mapping file
    id_col = mapping_df.columns[0]
    name_col = mapping_df.columns[1]
    
    # Clean up gene IDs (remove version numbers if present)
    mapping_df[id_col] = mapping_df[id_col].str.split('.').str[0]
    
    return dict(zip(mapping_df[id_col], mapping_df[name_col]))


def download_gtf_file(output_dir=None):
    """
    Download GENCODE v37 GTF file.
    
    Args:
        output_dir: Optional directory to save the file (uses current directory if not provided)
        
    Returns:
        Path to downloaded file
    """
    if output_dir:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        local_filename = os.path.join(output_dir, "gencode.v37.annotation.gtf.gz")
    else:
        local_filename = "gencode.v37.annotation.gtf.gz"
    
    print(f"Downloading GENCODE v37 GTF file from {GENCODE_V37_URL}")
    print(f"This file is large (~1GB) and may take some time to download...")
    
    # Simple progress indicator
    def report_progress(blocknum, blocksize, totalsize):
        percent = min(int(blocknum * blocksize * 100 / totalsize), 100)
        sys.stdout.write(f"\rDownload progress: {percent}% ({blocknum * blocksize / 1024 / 1024:.1f} MB)")
        sys.stdout.flush()
    
    try:
        urlretrieve(GENCODE_V37_URL, local_filename, reporthook=report_progress)
        print("\nDownload complete!")
        return local_filename
    except Exception as e:
        print(f"\nError downloading file: {str(e)}")
        sys.exit(1)


def extract_gene_id(attribute_str):
    """Extract gene ID from the attribute column of GTF file."""
    match = re.search(r'gene_id "([^"]+)"', attribute_str)
    if match:
        # Return gene ID without version if present
        gene_id = match.group(1)
        return gene_id.split('.')[0]
    return None


def extract_gene_name(attribute_str):
    """Extract gene name from the attribute column of GTF file."""
    match = re.search(r'gene_name "([^"]+)"', attribute_str)
    if match:
        return match.group(1)
    return None


def extract_gene_lengths(gtf_file, output_file):
    """
    Extract gene lengths from GENCODE GTF file.
    
    Args:
        gtf_file: Path to GENCODE GTF file (can be gzipped)
        output_file: Path to save the extracted gene lengths
        
    Returns:
        Path to the output file with gene lengths
    """
    print(f"Extracting gene lengths from {gtf_file}...")
    
    # Check if file is gzipped
    is_gzipped = gtf_file.endswith('.gz')
    
    # If gzipped, create a temporary unzipped file
    temp_file = None
    if is_gzipped:
        temp_file = gtf_file[:-3]  # Remove .gz extension
        print(f"Decompressing gzipped file to {temp_file}...")
        
        try:
            # Try using gunzip for better performance
            subprocess.run(['gunzip', '-c', gtf_file], stdout=open(temp_file, 'w'), check=True)
            gtf_file = temp_file
            print("Decompression complete!")
        except (subprocess.SubprocessError, FileNotFoundError):
            print("Could not use gunzip, falling back to Python's gzip module...")
            # Fallback to Python's gzip
            with gzip.open(gtf_file, 'rt') as f_in:
                with open(temp_file, 'w') as f_out:
                    for line in f_in:
                        f_out.write(line)
            gtf_file = temp_file
            print("Decompression complete!")
    
    # Process GTF file line by line
    gene_data = {}
    line_count = 0
    start_time = time.time()
    
    with open(gtf_file, 'r') as f:
        for line_num, line in enumerate(f):
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            line_count += 1
            
            # Print progress every million lines
            if line_num % 1000000 == 0 and line_num > 0:
                elapsed = time.time() - start_time
                rate = line_num / elapsed if elapsed > 0 else 0
                print(f"Processed {line_num:,} lines... ({rate:.1f} lines/sec)")
            
            # Parse line
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type != 'gene':
                continue
            
            # Extract gene info
            gene_id = extract_gene_id(fields[8])
            gene_name = extract_gene_name(fields[8])
            start = int(fields[3])
            end = int(fields[4])
            length = end - start + 1
            
            if gene_id:
                gene_data[gene_id] = {
                    'gene_name': gene_name,
                    'length': length
                }
    
    # Create DataFrame and save to TSV
    print(f"Creating gene length file: {output_file}")
    df = pd.DataFrame([
        {'gene_id': gene_id, 'gene_name': data['gene_name'], 'length': data['length']}
        for gene_id, data in gene_data.items()
    ])
    
    # Sort by gene_id
    df = df.sort_values('gene_id')
    
    # Save to file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Extracted lengths for {len(df):,} genes")
    
    # Print some statistics
    print(f"Length statistics:")
    print(f"  Min length: {df['length'].min():,} bp")
    print(f"  Max length: {df['length'].max():,} bp")
    print(f"  Median length: {df['length'].median():,.0f} bp")
    print(f"  Mean length: {df['length'].mean():,.0f} bp")
    
    # Clean up temporary file if created
    if temp_file and os.path.exists(temp_file) and is_gzipped:
        print(f"Removing temporary file: {temp_file}")
        try:
            os.remove(temp_file)
        except OSError as e:
            print(f"Warning: Could not remove temporary file: {str(e)}")
    
    return output_file


def load_gene_lengths(file_path):
    """
    Load gene length information from file.
    
    Args:
        file_path: Path to gene length file (tab-separated)
        
    Returns:
        Dictionary mapping gene IDs to gene lengths
    """
    print(f"Loading gene length information from {file_path}...")
    try:
        length_df = pd.read_csv(file_path, sep='\t')
        
        # Identify columns
        if len(length_df.columns) < 2:
            raise ValueError(f"Gene length file should have at least 2 columns (gene_id, length)")
        
        # Use the first column as gene_id and second as length if not specified
        id_col = 'gene_id' if 'gene_id' in length_df.columns else length_df.columns[0]
        length_col = 'length' if 'length' in length_df.columns else length_df.columns[1]
        
        # Clean up gene IDs (remove version numbers if present)
        if length_df[id_col].dtype == 'object':
            length_df[id_col] = length_df[id_col].str.split('.').str[0]
        
        # Ensure length column is numeric
        length_df[length_col] = pd.to_numeric(length_df[length_col], errors='coerce')
        
        # Create dictionary
        lengths = dict(zip(length_df[id_col], length_df[length_col]))
        print(f"Loaded length information for {len(lengths):,} genes")
        return lengths
    
    except Exception as e:
        print(f"Error loading gene lengths: {str(e)}")
        print("Using default gene length of 1000 bp for all genes")
        return {}  # Return empty dict, will use default length in normalization function


def normalize_to_tpm(rnaseq_df, gene_lengths, id_col='ID'):
    """
    Normalize raw counts to TPM (Transcripts Per Million) values.
    
    Args:
        rnaseq_df: DataFrame with RNA-seq count data
        gene_lengths: Dictionary or Series mapping gene IDs to gene lengths
        id_col: Column name for gene IDs
        
    Returns:
        DataFrame with normalized TPM values
    """
    print("Normalizing raw counts to TPM values...")
    
    # Make a copy to avoid modifying the original
    tpm_df = rnaseq_df.copy()
    
    # Get count columns (all columns except ID)
    count_cols = [col for col in tpm_df.columns if col != id_col]
    
    # Extract gene IDs and ensure they're in the same format as the length dictionary
    gene_ids = tpm_df[id_col]
    if gene_ids.dtype == 'object' and isinstance(gene_ids.iloc[0], str):
        # Remove version numbers if present
        gene_ids_no_version = gene_ids.str.split('.').str[0]
    else:
        gene_ids_no_version = gene_ids
    
    # Create a Series of gene lengths for each gene in our data
    lengths = pd.Series([gene_lengths.get(gene_id, 1000) for gene_id in gene_ids_no_version], 
                       index=tpm_df.index)
    
    # Print warning for genes with missing lengths
    missing_lengths = sum(1 for gene_id in gene_ids_no_version if gene_id not in gene_lengths)
    if missing_lengths > 0:
        print(f"Warning: {missing_lengths:,} genes are missing length information. Using default length of 1000 bp.")
    
    # Calculate TPM for each column
    for col in count_cols:
        # Step 1: Divide counts by gene length to get reads per kilobase (RPK)
        rpk = tpm_df[col] / (lengths / 1000)
        
        # Step 2: Calculate scaling factor (sum of RPK / 1,000,000)
        scaling_factor = rpk.sum() / 1e6
        
        # Step 3: Divide RPK by scaling factor to get TPM
        if scaling_factor > 0:  # Avoid division by zero
            tpm_df[col] = rpk / scaling_factor
        else:
            print(f"Warning: Zero scaling factor for {col}. Setting all TPM values to 0.")
            tpm_df[col] = 0
    
    # Print some statistics
    for col in count_cols:
        print(f"  {col}: Min TPM = {tpm_df[col].min():.2f}, Max TPM = {tpm_df[col].max():.2f}, Mean TPM = {tpm_df[col].mean():.2f}")
    
    return tpm_df


def assign_quantiles(data, condition, min_expr=1.0, num_quantiles=10):
    """
    Assign genes to quantiles based on expression in a specific condition.
    
    Args:
        data: DataFrame with RNA-seq data
        condition: Column name for the condition to analyze
        min_expr: Minimum expression threshold
        num_quantiles: Number of quantiles to create
        
    Returns:
        Series with gene IDs as index and quantile groups as values
    """
    # Create a series of expression values
    expr = data[condition].copy()
    
    # For quantile assignment, treat non-expressed genes specially
    non_expressed = expr < min_expr
    
    # Assign quantile 0 to non-expressed genes (changed from 1)
    quantiles = pd.Series(0, index=expr.index)
    
    # For expressed genes, create num_quantiles quantiles (1 through num_quantiles)
    if (~non_expressed).sum() > 0:
        expr_genes = expr[~non_expressed]
        # Use qcut for equal-sized quantiles based on expression values
        try:
            # Create n quantiles (labels 1-10 instead of 2-10)
            n_quantiles = num_quantiles  # Number of quantiles for expressed genes
            q_labels = list(range(1, num_quantiles + 1))  # Changed from 2-11 to 1-10
            
            # Generate quantiles with proper number of bins
            expr_quantiles = pd.qcut(expr_genes, 
                                    n_quantiles, 
                                    labels=q_labels,
                                    duplicates='drop')  # Handle duplicate edges
            
            # Update quantile assignments for expressed genes
            quantiles.loc[expr_quantiles.index] = expr_quantiles
        except ValueError as e:
            print(f"Warning for {condition}: {str(e)}")
            print(f"Trying alternative quantile approach...")
            
            # Alternative approach using percentiles
            try:
                # Calculate rank percentiles
                ranks = expr_genes.rank(pct=True)
                # Scale to our quantile range (1 to num_quantiles)
                scaled_ranks = 1 + (ranks * (num_quantiles - 1)).round().astype(int)  # Changed from 2-10 to 1-10
                # Ensure no values exceed the max quantile
                scaled_ranks = scaled_ranks.clip(upper=num_quantiles)
                # Update quantiles
                quantiles.loc[scaled_ranks.index] = scaled_ranks
            except Exception as e2:
                print(f"Alternative approach failed: {str(e2)}")
                # Last resort: simple binning by value
                min_val = expr_genes.min()
                max_val = expr_genes.max()
                if min_val == max_val:
                    # All expressed genes have same value, assign to highest quantile
                    quantiles.loc[expr_genes.index] = num_quantiles
                else:
                    # Create bins by value range
                    bin_edges = np.linspace(min_val, max_val, num_quantiles)
                    for i, gene_id in enumerate(expr_genes.index):
                        # Find which bin this gene belongs to
                        bin_idx = np.digitize(expr_genes.iloc[i], bin_edges)
                        # Convert to our quantile scale (1 to num_quantiles)
                        q_val = min(bin_idx, num_quantiles)  # Changed to start from 1
                        quantiles.loc[gene_id] = q_val
    
    return quantiles


def assign_equal_bins_quantiles(data, condition, min_expr=1.0, num_quantiles=10):
    """
    Assign genes to strictly equal-sized quantiles for expressed genes only.
    
    Args:
        data: DataFrame with RNA-seq data
        condition: Column name for the condition to analyze
        min_expr: Minimum expression threshold
        num_quantiles: Number of quantiles to create
        
    Returns:
        Series with gene IDs as index and quantile groups as values
    """
    # Create a series of expression values
    expr = data[condition].copy()
    
    # Separate non-expressed and expressed genes
    non_expressed = expr < min_expr
    expressed = ~non_expressed
    
    # Start with all genes assigned to quantile 0 (non-expressed) - changed from 1
    quantiles = pd.Series(0, index=expr.index)
    
    # Process expressed genes only if there are any
    if expressed.sum() > 0:
        expr_genes = expr[expressed]
        
        # Sort the expressed genes by expression value
        sorted_genes = expr_genes.sort_values()
        
        # Calculate how many genes should be in each quantile
        n_expressed = len(sorted_genes)
        genes_per_quantile = n_expressed / num_quantiles  # Changed from num_quantiles-1
        
        # Assign quantiles by equal gene counts rather than by expression value
        for q in range(1, num_quantiles + 1):  # Changed from 2-11 to 1-10
            # Calculate the start and end index for this quantile
            start_idx = int((q - 1) * genes_per_quantile)  # Changed to index from 0
            end_idx = int(q * genes_per_quantile)
            # For the last quantile, ensure we include all remaining genes
            if q == num_quantiles:
                end_idx = n_expressed
            
            # Get the gene IDs for this quantile
            if start_idx < end_idx:  # Ensure we have at least one gene
                quantile_genes = sorted_genes.index[start_idx:end_idx]
                # Assign the quantile number to these genes
                quantiles.loc[quantile_genes] = q
    
    return quantiles


def create_output_file(gene_ids, gene_names, quantiles, output_path):
    """
    Create output file with gene IDs, names, and quantile groups.
    
    Args:
        gene_ids: Series or list of gene IDs
        gene_names: Dictionary mapping gene IDs to gene names
        quantiles: Series with quantile group assignments
        output_path: Path to save output file
    """
    # Convert quantiles to integers - important for proper output format
    quantiles = quantiles.astype(int)
    
    # Create DataFrame with correct columns
    output_df = pd.DataFrame({
        'gene_id_no_version': gene_ids.values if hasattr(gene_ids, 'values') else gene_ids,
        'Group_Num': quantiles.values
    })
    
    # Add gene names if available
    if gene_names:
        output_df['gene_name'] = output_df['gene_id_no_version'].map(
            lambda x: gene_names.get(x, x)
        )
        # Reorder columns to match requested format
        output_df = output_df[['gene_id_no_version', 'gene_name', 'Group_Num']]
    
    # Save to file with tab separator and ensure proper formatting
    output_df.to_csv(output_path, sep='\t', index=False)
    
    # Verify file was created correctly
    try:
        # Read back the first few lines to check format
        with open(output_path, 'r') as f:
            first_lines = [next(f) for _ in range(min(5, len(output_df) + 1))]
        print(f"File format verification (first few lines):")
        for line in first_lines:
            print(f"  {line.strip()}")
    except Exception as e:
        print(f"Warning: Could not verify file format: {str(e)}")
    
    # Print distribution of quantiles
    quantile_counts = output_df['Group_Num'].value_counts().sort_index()
    print(f"Quantile distribution for {os.path.basename(output_path)}:")
    for group, count in quantile_counts.items():
        print(f"  Group {group}: {count:,} genes")


def process_rnaseq_data(args):
    """
    Main function to process RNA-seq data and generate quantile files.
    
    Args:
        args: Command line arguments
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"Created output directory: {args.output_dir}")
    
    # Create normalized data subdirectory
    norm_dir = os.path.join(args.output_dir, "normalized_data")
    if not os.path.exists(norm_dir):
        os.makedirs(norm_dir)
        print(f"Created normalized data directory: {norm_dir}")
    
    # Create equal bins subdirectory if needed and not in matrix-only mode
    equal_bins_dir = None
    if args.equal_bins and not args.normalize_matrix_only:
        equal_bins_dir = os.path.join(args.output_dir, "equal_bins")
        if not os.path.exists(equal_bins_dir):
            os.makedirs(equal_bins_dir)
            print(f"Created equal bins directory: {equal_bins_dir}")
    
    # Create gene length directory if normalization is requested
    gene_length_dir = None
    gene_length_file = None
    if args.normalize or args.normalize_matrix_only or not args.gene_name_file:
        gene_length_dir = os.path.join(args.output_dir, "gene_length_data")
        if not os.path.exists(gene_length_dir):
            os.makedirs(gene_length_dir)
            print(f"Created gene length directory: {gene_length_dir}")
    
    # Load RNA-seq data
    rnaseq_df = load_rnaseq_data(args.input, args.id_col)
    
    # Process gene length and name data
    gene_lengths = {}
    gene_names = {}
    
    # Determine if we need to extract gene information from GENCODE
    need_gencode_data = args.normalize or args.normalize_matrix_only or not args.gene_name_file
    
    if need_gencode_data:
        # Check if gene length file is provided
        if args.gene_length_file and os.path.exists(args.gene_length_file):
            # Use provided gene length file
            gene_length_file = args.gene_length_file
            gene_lengths = load_gene_lengths(gene_length_file)
            
            # Also extract gene names from the same file if it has them
            try:
                length_df = pd.read_csv(gene_length_file, sep='\t')
                if 'gene_id' in length_df.columns and 'gene_name' in length_df.columns:
                    gene_names = dict(zip(
                        length_df['gene_id'].str.split('.').str[0], 
                        length_df['gene_name']
                    ))
                    print(f"Extracted {len(gene_names):,} gene names from gene length file")
            except Exception as e:
                print(f"Could not extract gene names from gene length file: {str(e)}")
        else:
            # Need to extract gene lengths from GENCODE
            gtf_file = args.gtf_file
            
            if not gtf_file or not os.path.exists(gtf_file):
                if args.no_download:
                    print("ERROR: Gene information needed but no gene length file or GTF file provided,")
                    print("and --no_download flag is set. Cannot proceed.")
                    sys.exit(1)
                
                # Download GTF file
                gtf_file = download_gtf_file(gene_length_dir)
            
            # Extract gene lengths from GTF file
            gene_length_file = os.path.join(gene_length_dir, "gencode_v37_gene_lengths.tsv")
            extract_gene_lengths(gtf_file, gene_length_file)
            
            # Load extracted gene lengths
            gene_lengths = load_gene_lengths(gene_length_file)
            
            # Extract gene names from the same file
            try:
                length_df = pd.read_csv(gene_length_file, sep='\t')
                gene_names = dict(zip(
                    length_df['gene_id'].str.split('.').str[0], 
                    length_df['gene_name']
                ))
                print(f"Extracted {len(gene_names):,} gene names from GENCODE data")
            except Exception as e:
                print(f"Could not extract gene names from GENCODE data: {str(e)}")
    
    # If gene_name_file is specified, use it (overrides names from gene length file)
    if args.gene_name_file and os.path.exists(args.gene_name_file):
        gene_names = load_gene_names(args.gene_name_file)
    
    # Normalize to TPM if requested or in matrix-only mode
    if args.normalize or args.normalize_matrix_only:
        # Normalize counts to TPM
        tpm_df = normalize_to_tpm(rnaseq_df, gene_lengths, args.id_col)
        print("Normalized raw counts to TPM values")
        
        # Save normalized matrix (full dataframe)
        norm_matrix_file = os.path.join(norm_dir, "normalized_tpm_matrix.tsv")
        tpm_df.to_csv(norm_matrix_file, sep='\t', index=False)
        print(f"Saved full normalized TPM matrix to {norm_matrix_file}")
        
        # If matrix-only mode, we're done
        if args.normalize_matrix_only:
            print("\nTPM normalization completed successfully. Matrix saved to:")
            print(f"  {norm_matrix_file}")
            return
        
        # Otherwise continue with the regular workflow using normalized data
        rnaseq_df = tpm_df
    
    # Extract gene IDs without version numbers if they have them
    gene_ids = rnaseq_df[args.id_col].copy()
    if gene_ids.dtype == 'object':
        gene_ids = gene_ids.str.split('.').str[0]
    
    # Identify and process each condition (all columns except the ID column)
    conditions = [col for col in rnaseq_df.columns if col != args.id_col]
    
    for condition in conditions:
        print(f"\nProcessing condition: {condition}")
        
        # Set index for calculation but keep a copy of the original dataframe for output
        indexed_df = rnaseq_df.set_index(args.id_col)
        
        # Calculate standard quantiles for this condition
        quantiles = assign_quantiles(
            indexed_df, 
            condition, 
            args.min_expr, 
            args.quantiles
        )
        
        # Create output file name based on condition
        output_file = os.path.join(
            args.output_dir, 
            f"{condition}_quantiles.tsv"
        )
        
        # Create and save output file
        create_output_file(gene_ids, gene_names, quantiles, output_file)
        
        # If equal bins option is enabled, create those files too
        if args.equal_bins:
            print(f"Creating equal-sized bins for expressed genes in {condition}")
            
            # Calculate equal-sized quantiles
            equal_quantiles = assign_equal_bins_quantiles(
                indexed_df,
                condition,
                args.min_expr,
                args.quantiles
            )
            
            # Create output file for equal bins
            equal_output_file = os.path.join(
                equal_bins_dir,
                f"{condition}_equal_quantiles.tsv"
            )
            
            # Create and save output file
            create_output_file(gene_ids, gene_names, equal_quantiles, equal_output_file)


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Enable normalization if normalize_matrix_only is set
    if args.normalize_matrix_only:
        args.normalize = True
    
    # Validate arguments
    if args.normalize and not args.gene_length_file and not args.gtf_file and args.no_download:
        print("ERROR: Normalization requires gene length information.")
        print("Please provide either --gene_length_file, --gtf_file, or allow automatic download.")
        sys.exit(1)
    
    process_rnaseq_data(args)
    
    if not args.normalize_matrix_only:
        print("\nRNA-seq quantile analysis completed successfully!")
        if args.normalize:
            print("Data was normalized to TPM values before quantile assignment.")
            print("Full normalized matrix was saved to normalized_data/normalized_tpm_matrix.tsv")
        else:
            print("Raw counts were used for quantile assignment.")
        
        print("\nOutput files contain:")
        print("- gene_id_no_version: Ensembl gene ID without version number")
        if args.gene_name_file or args.normalize or args.gtf_file:
            print("- gene_name: Gene symbol extracted from GENCODE or provided gene name file")
        print("- Group_Num: Quantile group (0 for non-expressed, 1-10 for expressed genes)")
        print("\nTo view the results, you can use commands like:")
        print(f"  head {os.path.join(args.output_dir, '[condition]_quantiles.tsv')}")


if __name__ == "__main__":
    main()
