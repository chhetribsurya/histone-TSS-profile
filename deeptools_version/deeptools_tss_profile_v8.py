#!/usr/bin/env python3
import pandas as pd
import re
import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def download_gencode_gtf(version="37"):
    """Download Gencode GTF file if not already present"""
    gtf_file = f"gencode.v{version}.annotation.gtf"
    
    if not os.path.exists(gtf_file) and not os.path.exists(f"{gtf_file}.gz"):
        print(f"Downloading {gtf_file}...")
        download_cmd = f"wget -O {gtf_file}.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/{gtf_file}.gz"
        subprocess.run(download_cmd, shell=True, check=True)
    
    # Decompress if needed
    if not os.path.exists(gtf_file) and os.path.exists(f"{gtf_file}.gz"):
        print(f"Decompressing {gtf_file}.gz...")
        decompress_cmd = f"gunzip -k {gtf_file}.gz"
        subprocess.run(decompress_cmd, shell=True, check=True)
    
    return gtf_file

def download_chrom_sizes(genome="hg37"):
    """Download chromosome sizes file if not already present"""
    chrom_sizes_file = f"{genome}.chrom.sizes"
    
    if not os.path.exists(chrom_sizes_file):
        print(f"Downloading {chrom_sizes_file}...")
        download_cmd = f"wget -O {chrom_sizes_file} https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/{chrom_sizes_file}"
        subprocess.run(download_cmd, shell=True, check=True)
    
    return chrom_sizes_file

def convert_wig_to_bigwig(wig_file, chrom_sizes_file):
    """Convert WIG file to BigWig format"""
    bigwig_file = wig_file
    if wig_file.endswith(".gz"):
        bigwig_file = wig_file[:-3]  # Remove .gz
    
    bigwig_file = bigwig_file.replace(".wig", ".bw")
    
    if not os.path.exists(bigwig_file):
        print(f"Converting {wig_file} to {bigwig_file}...")
        
        # If file is compressed, decompress first
        if wig_file.endswith(".gz"):
            uncompressed_wig = wig_file[:-3]
            if not os.path.exists(uncompressed_wig):
                print(f"Decompressing {wig_file}...")
                decompress_cmd = f"gunzip -k {wig_file}"
                subprocess.run(decompress_cmd, shell=True, check=True)
            wig_file = uncompressed_wig
        
        # Convert WIG to BigWig
        convert_cmd = f"wigToBigWig {wig_file} {chrom_sizes_file} {bigwig_file}"
        subprocess.run(convert_cmd, shell=True, check=True)
    else:
        print(f"Bigwig file exists: {bigwig_file}...")
    
    return bigwig_file

def create_tss_bed_file(genes_df, gtf_file, output_dir="results", output_bed="gene_tss_with_groups.bed"):
    """Create BED file with TSS positions and groups"""
    print("Creating TSS BED file with groups...")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Full path for output bed file
    output_bed_path = os.path.join(output_dir, output_bed)
    
    # Create a mapping from gene ID to group
    gene_groups = dict(zip(genes_df.gene_id_no_version, genes_df.Group_Num))
    
    # Create separate BED files for each group
    group_beds = {}
    
    # Parse GTF to extract TSS positions
    tss_positions = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "gene":
                continue
                
            # Extract gene_id without version
            gene_id_match = re.search(r'gene_id "(ENSG\d+)\.\d+"', fields[8])
            if not gene_id_match:
                continue
                
            gene_id = gene_id_match.group(1)
            
            # Check if this gene is in our list
            if gene_id not in gene_groups:
                continue
                
            # Get gene name if available
            gene_name_match = re.search(r'gene_name "([^"]+)"', fields[8])
            gene_name = gene_name_match.group(1) if gene_name_match else gene_id
                
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            group = gene_groups[gene_id]
            
            # TSS is at start for + strand genes, end for - strand genes
            tss = start if strand == "+" else end
            
            tss_positions.append([chrom, tss-1, tss, f"{gene_id}_{gene_name}", 0, strand, group])
    
    # Create unique groups
    unique_groups = set([pos[6] for pos in tss_positions])
    
    # Write master BED file with group in 6th column
    with open(output_bed_path, "w") as f:
        for pos in tss_positions:
            f.write("\t".join(map(str, pos)) + "\n")
    
    # Create separate BED files for each group
    for group in unique_groups:
        group_file = os.path.join(output_dir, f"group_{group}_tss.bed")
        with open(group_file, "w") as f:
            for pos in tss_positions:
                if pos[6] == group:
                    # Write without the group column
                    f.write("\t".join(map(str, pos[:6])) + "\n")
    
    # Create "All" group BED file
    all_file = os.path.join(output_dir, "group_All_tss.bed")
    with open(all_file, "w") as f:
        for pos in tss_positions:
            f.write("\t".join(map(str, pos[:6])) + "\n")
    
    # Create a sample file that lists all group BED files
    sample_file = os.path.join(output_dir, "group_samples.txt")
    with open(sample_file, "w") as f:
        for group in unique_groups:
            f.write(f"{output_dir}/group_{group}_tss.bed\t{group}\n")
        f.write(f"{all_file}\tAll\n")
    
    print(f"Created BED file: {output_bed_path}")
    print(f"Created sample file: {sample_file}")
    return sample_file

def run_deeptools(bigwig_file, sample_file, output_dir="results", output_prefix="H1.3"):
    """Run deepTools computeMatrix and plotProfile using multiple BED files"""
    print("Running deepTools computeMatrix...")
    
    # Run computeMatrix in multiple mode
    matrix_file = os.path.join(output_dir, f"{output_prefix}_matrix.gz")
    compute_cmd = (
        f"computeMatrix reference-point "
        f"--referencePoint TSS "
        f"--scoreFileName {bigwig_file} "
        f"--regionsFileName {output_dir}/group_*_tss.bed "
        f"--beforeRegionStartLength 3000 "
        f"--afterRegionStartLength 3000 "
        f"--binSize 1 "
        f"--sortRegions keep "
        f"--smartLabels "
        f"-o {matrix_file}"
    )
    subprocess.run(compute_cmd, shell=True, check=True)
    
    # Run plotProfile
    plot_file = os.path.join(output_dir, f"{output_prefix}_TSS_profile.pdf")
    #f"--colors 'grey' '#FF0000' '#FF4500' '#FFA500' '#FFD700' '#ADFF2F' '#32CD32' '#20B2AA' '#1E90FF' '#4169E1' '#8A2BE2' '#000000' "
    profile_cmd = (
        f"plotProfile "
        f"--matrixFile {matrix_file} "
        f"--outFileName {plot_file} "
        f"--plotTitle \"{output_prefix}\" "
        f"--refPointLabel \"TSS\" "
        f"--yAxisLabel \"Average Profile\" "
        f"--colors '#FF0000' '#FF4500' '#FFA500' '#FFD700' '#ADFF2F' '#32CD32' '#20B2AA' '#1E90FF' '#4169E1' '#8A2BE2' '#000000' "
        f"--legendLocation \"lower-left\" "
        f"--plotHeight 16 "
        f"--plotWidth 12 "
    )
    subprocess.run(profile_cmd, shell=True, check=True)
    
    print(f"Created plot file: {plot_file}")

def main():
    # Check if required tools are installed
    try:
        subprocess.run("which wigToBigWig", shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error: wigToBigWig not found. Please install UCSC tools.")
        sys.exit(1)
    
    try:
        subprocess.run("which computeMatrix", shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error: computeMatrix not found. Please install deepTools.")
        sys.exit(1)
    
    # Process command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Generate TSS profile plots')
    parser.add_argument('--wig', type=str, required=True, help='WIG file path')
    parser.add_argument('--genes', type=str, required=True, help='Gene list TSV file path')
    parser.add_argument('--gtf_version', type=str, default='37', help='Gencode GTF version (default: 37)')
    parser.add_argument('--genome', type=str, default='hg37', help='Genome build (default: hg37)')
    parser.add_argument('--output_prefix', type=str, default='H1.3', help='Output prefix for files')
    parser.add_argument('--groups', type=str, default='1,2,3,4,5,All', 
                    help='Comma-separated Group_Num values to include (default: "1,2,3,4,5,All")')
    parser.add_argument('--output_dir', type=str, default='results', 
                    help='Directory to save output files (default: "results")')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Download necessary files
    gtf_file = download_gencode_gtf(args.gtf_version)
    chrom_sizes_file = download_chrom_sizes(args.genome)
    
    # Read gene list
    print(f"Reading gene list from {args.genes}...")
    try:
        genes_df = pd.read_csv(args.genes, sep='\t')

        # Process selected groups properly (handle numeric and string 'All')
        selected_groups_input = [g.strip() for g in args.groups.split(',')]
        numeric_groups = [int(g) for g in selected_groups_input if g.isdigit()]
        string_groups = [g for g in selected_groups_input if not g.isdigit()]

        # subset dataframe properly based on Group_Num types
        subset_numeric = genes_df[genes_df['Group_Num'].isin(numeric_groups)] if numeric_groups else pd.DataFrame()
        subset_string = genes_df[genes_df['Group_Num'].astype(str).isin(string_groups)] if string_groups else pd.DataFrame()
        # Combine subsets
        genes_df = pd.concat([subset_numeric, subset_string])

        print(f"Subset to groups: {selected_groups_input}")
        print(genes_df)
        
        # Check that required columns exist
        required_columns = ['gene_id_no_version', 'Group_Num']
        missing_columns = [col for col in required_columns if col not in genes_df.columns]
        
        if missing_columns:
            print(f"Error: Missing required columns in gene list file: {', '.join(missing_columns)}")
            print(f"Available columns: {', '.join(genes_df.columns)}")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error reading gene list file: {e}")
        sys.exit(1)
    
    # Convert WIG to BigWig
    try:
        bigwig_file = convert_wig_to_bigwig(args.wig, chrom_sizes_file)
    except Exception as e:
        print(f"Error converting WIG to BigWig: {e}")
        sys.exit(1)
    
    # Create TSS BED file
    try:
        bed_file = create_tss_bed_file(genes_df, gtf_file, args.output_dir)
    except Exception as e:
        print(f"Error creating TSS BED file: {e}")
        sys.exit(1)
    
    # Run deepTools
    try:
        run_deeptools(bigwig_file, bed_file, args.output_dir, args.output_prefix)
    except Exception as e:
        print(f"Error running deepTools: {e}")
        sys.exit(1)
    
    print("Done!")

if __name__ == "__main__":
    main()
