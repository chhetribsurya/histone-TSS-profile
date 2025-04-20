#!/usr/bin/env python3
"""
Script to create TSS BED files from RefSeq gene annotations and expression group files.
This script:
1. Downloads RefSeq gene annotations if not available
2. Reads expression group files containing RefSeq IDs
3. Generates a TSS BED file with expression groups
4. Outputs individual BED files for each expression group

Advanced features include:
- Multiple TSS selection methods for handling alternative transcripts
- Support for different genome builds
- Detailed statistics reporting
"""

import os
import sys
import argparse
import re
import subprocess
import gzip
from collections import defaultdict

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Generate TSS BED files from RefSeq annotations and expression groups')
    parser.add_argument('--refseq', type=str, help='Path to RefSeq gene annotation file (will download if not provided)')
    parser.add_argument('--expr_dir', type=str, required=True, help='Directory containing expression group files')
    parser.add_argument('--output_dir', type=str, default='results', help='Directory to save output BED files')
    parser.add_argument('--output_bed', type=str, default='refseq_tss_with_groups.bed', help='Filename for the main output BED file')
    parser.add_argument('--genome', type=str, default='hg38', help='Genome build (default: hg38)')
    parser.add_argument('--verbose', action='store_true', help='Print detailed debug information')
    parser.add_argument('--tss_selection', type=str, 
                      choices=['canonical', 'most_upstream', 'most_downstream', 'longest_transcript', 
                               'shortest_transcript', 'first', 'all'], 
                      default='canonical', 
                      help='Method to select representative TSS: canonical (default), most_upstream, most_downstream, '
                           'longest_transcript, shortest_transcript, first, or all (retain all TSS positions)')
    return parser.parse_args()

def download_refseq_genes(output_dir, genome="hg38"):
    """Download RefSeq gene annotation file if not already present, save in a subdirectory"""
    # Create refseq_genefiles subdirectory
    refseq_dir = os.path.join(output_dir, "refseq_genefiles")
    os.makedirs(refseq_dir, exist_ok=True)
    
    refseq_file = os.path.join(refseq_dir, f"refGene_{genome}.txt.gz")
    
    if not os.path.exists(refseq_file):
        print(f"Downloading RefSeq gene annotations for {genome}...")
        download_cmd = f"wget -O {refseq_file} https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/refGene.txt.gz"
        try:
            subprocess.run(download_cmd, shell=True, check=True)
            print(f"Downloaded RefSeq annotations to {refseq_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading RefSeq annotations: {e}")
            sys.exit(1)
    else:
        print(f"Using existing RefSeq annotations: {refseq_file}")
    
    return refseq_file

def read_expression_group_files(expr_dir, verbose=False):
    """
    Read expression group files and create a mapping from RefSeq ID to group.
    Each file should be named with its group number preceded by "NM" (e.g., "T47D_RNAseq_GeneName_NM_0.txt", i.e. "*NM_0.txt").
    Files contain RefSeq IDs (NM_*, NR_*), one per line.
    """
    print(f"Reading expression group files from {expr_dir}...")
    
    gene_to_group = {}
    group_stats = defaultdict(int)
    groups_found = []
    total_refseq_count = 0
    
    # Process files in a specific order to ensure deterministic behavior
    filenames = sorted(os.listdir(expr_dir))
    
    for filename in filenames:
        filepath = os.path.join(expr_dir, filename)
        if not os.path.isfile(filepath):
            continue
        
        # Extract group number from filename
        group = None
        if filename == "All.txt":
            # Skip the "All.txt" file as requested
            print(f"Skipping All.txt file as requested")
            continue
        else:
            match = re.search(r'NM_(\d+)\.txt$', filename)
            if match:
                group = int(match.group(1))  # Make sure group is stored as int
        
        if group is None:
            continue
            
        groups_found.append(group)
        count = 0
        
        with open(filepath, 'r') as f:
            for line in f:
                refseq_id = line.strip()
                if refseq_id:
                    # Store multiple groups per RefSeq ID if needed (handle duplicates)
                    gene_to_group[refseq_id] = group
                    group_stats[group] += 1
                    count += 1
                    total_refseq_count += 1
        
        print(f"Reading group {group} file: {filename} - {count} RefSeq IDs")
    
    # Count unique RefSeq IDs
    unique_refseq_ids = len(gene_to_group)
    
    print(f"\nFound {len(groups_found)} expression groups: {sorted(groups_found, key=lambda x: (0 if isinstance(x, str) else 1, x if isinstance(x, str) else x))}")
    print(f"Read {total_refseq_count} total RefSeq IDs from files")
    print(f"Unique RefSeq IDs: {unique_refseq_ids}")
    
    if verbose:
        print("\nRefSeq IDs per group:")
        for group in sorted(group_stats.keys(), key=lambda x: (0 if isinstance(x, str) else 1, x if isinstance(x, str) else x)):
            print(f"  Group {group}: {group_stats[group]} RefSeq IDs")
    
    # Return both the mapping and the original count statistics for each group
    return gene_to_group, group_stats

def parse_refseq_genes(refseq_file, gene_to_group, original_group_counts, args):
    """
    Parse RefSeq gene annotations and extract TSS positions for genes in the expression groups.
    Apply the selected method to handle multiple transcripts for the same RefSeq ID.
    RefSeq format: bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, etc.
    """
    print(f"Parsing RefSeq gene annotations from {refseq_file}...")
    verbose = args.verbose
    
    # List to store TSS positions
    tss_positions = []
    refseq_ids_found = set()
    group_counts = defaultdict(int)
    
    # Count total RefSeq entries in the annotation file
    total_refseq_entries = 0
    
    # Track duplicate RefSeq entries
    refseq_duplicates = defaultdict(int)
    
    # Store complete entries for selection
    refseq_entries = {}
    
    # Open file (handle gzipped files)
    opener = gzip.open if refseq_file.endswith('.gz') else open
    
    # First pass - count entries and duplicates
    with opener(refseq_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 16:  # Ensure we have all necessary fields
                name = fields[1]
                # Only count entries with NM_ or NR_ prefixes
                if name.startswith("NM_") or name.startswith("NR_"):
                    total_refseq_entries += 1
                    refseq_duplicates[name] += 1
    
    print(f"Total RefSeq entries in annotation file: {total_refseq_entries}")
    
    # Check for duplicates in the RefSeq file
    duplicate_count = sum(1 for name, count in refseq_duplicates.items() if count > 1)
    if duplicate_count > 0:
        print(f"Found {duplicate_count} RefSeq IDs with multiple entries in the annotation file")
        if verbose:
            for name, count in sorted(refseq_duplicates.items(), key=lambda x: x[1], reverse=True):
                if count > 1 and count < 5:  # Limit output for readability
                    print(f"  RefSeq ID {name} appears {count} times in the annotation file")
    
    # Second pass - parse and store all entries
    with opener(refseq_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            
            # Check if it's a valid RefSeq line
            if len(fields) < 16:  # Need cdsStart, cdsEnd, exonCount, etc.
                continue
                
            try:
                name = fields[1]  # RefSeq ID (NM_* or NR_*)
                chrom = fields[2]
                strand = fields[3]
                tx_start = int(fields[4])
                tx_end = int(fields[5])
                cds_start = int(fields[6])
                cds_end = int(fields[7])
                exon_count = int(fields[8])
                exon_starts = fields[9].strip(',').split(',')
                exon_ends = fields[10].strip(',').split(',')
                gene_symbol = fields[12]  # name2 field is the gene symbol
                
                # Skip non-standard chromosomes
                if '_' in chrom or chrom.startswith('chrUn'):
                    continue
                
                # Check if this RefSeq ID is in our expression groups
                if name in gene_to_group:
                    # TSS is at txStart for + strand genes, txEnd for - strand genes
                    tss = tx_start if strand == "+" else tx_end
                    
                    # Store the group EXACTLY as it appears in the dictionary
                    # This is critical for maintaining type (int vs. str)
                    group = gene_to_group[name]
                    
                    # Create a key for the refseq ID + group combo
                    key = (name, group)
                    
                    # Calculate transcript length
                    transcript_length = tx_end - tx_start
                    
                    # Create a complete entry for selection
                    entry = {
                        'refseq_id': name,
                        'group': group,
                        'chrom': chrom,
                        'strand': strand,
                        'tx_start': tx_start,
                        'tx_end': tx_end,
                        'tss': tss,
                        'cds_start': cds_start,
                        'cds_end': cds_end,
                        'exon_count': exon_count,
                        'transcript_length': transcript_length,
                        'gene_symbol': gene_symbol,
                        'is_canonical': cds_start < cds_end and cds_start > 0 and cds_end > 0,  # Has valid CDS
                        'tss_position': [
                            chrom, tss-1, tss, f"{name}_{gene_symbol}", 0, strand, group
                        ]
                    }
                    
                    # Add to entries
                    if key not in refseq_entries:
                        refseq_entries[key] = []
                    refseq_entries[key].append(entry)
                    
                    # Count this as found
                    refseq_ids_found.add(name)
                    
            except Exception as e:
                if verbose:
                    print(f"Error parsing line: {e}")
                continue
    
    # Process entries based on command line args
    if args.tss_selection == 'all':
        print("\nRetaining all TSS positions as requested ('all' selection method)...")
        # Add all TSS positions
        for key, entries in refseq_entries.items():
            name, group = key
            for entry in entries:
                tss_positions.append(entry['tss_position'])
                group_counts[group] += 1
            if verbose and len(entries) > 1:
                print(f"  Retained all {len(entries)} TSS positions for {name} (Group {group})")
    else:
        # Select representative TSS positions
        print(f"\nSelecting representative TSS positions using '{args.tss_selection}' method...")
        
        for key, entries in refseq_entries.items():
            name, group = key
            
            # If there's only one entry, use it
            if len(entries) == 1:
                tss_positions.append(entries[0]['tss_position'])
                group_counts[group] += 1
                continue
            
            # For multiple entries, apply selection strategy
            selected_entry = None
            
            if args.tss_selection == 'canonical':
                # Prefer entries with valid CDS
                canonical_entries = [e for e in entries if e['is_canonical']]
                if canonical_entries:
                    # Among canonical, prefer the longest transcript
                    selected_entry = max(canonical_entries, key=lambda e: e['transcript_length'])
                else:
                    # If no canonical, use the longest transcript
                    selected_entry = max(entries, key=lambda e: e['transcript_length'])
            
            elif args.tss_selection == 'most_upstream':
                # Select the most upstream TSS relative to gene direction
                if entries[0]['strand'] == '+':
                    selected_entry = min(entries, key=lambda e: e['tss'])
                else:
                    selected_entry = max(entries, key=lambda e: e['tss'])
            
            elif args.tss_selection == 'most_downstream':
                # Select the most downstream TSS relative to gene direction
                if entries[0]['strand'] == '+':
                    selected_entry = max(entries, key=lambda e: e['tss'])
                else:
                    selected_entry = min(entries, key=lambda e: e['tss'])
            
            elif args.tss_selection == 'longest_transcript':
                # Select the entry with the longest transcript
                selected_entry = max(entries, key=lambda e: e['transcript_length'])
            
            elif args.tss_selection == 'shortest_transcript':
                # Select the entry with the shortest transcript
                selected_entry = min(entries, key=lambda e: e['transcript_length'])
            
            else:  # 'first' or default
                # Just take the first entry
                selected_entry = entries[0]
            
            # Add the selected entry
            tss_positions.append(selected_entry['tss_position'])
            group_counts[group] += 1
            
            if verbose and len(entries) > 1:
                print(f"  Selected 1 of {len(entries)} TSS positions for {name} (Group {group})")
    
    # Calculate percentage of RefSeq IDs found
    found_percentage = (len(refseq_ids_found) / len(gene_to_group)) * 100 if len(gene_to_group) > 0 else 0
    missing_count = len(gene_to_group) - len(refseq_ids_found)
    missing_percentage = (missing_count / len(gene_to_group)) * 100 if len(gene_to_group) > 0 else 0
    
    print(f"Found {len(refseq_ids_found)} out of {len(gene_to_group)} RefSeq IDs in the annotations ({found_percentage:.2f}%)")
    print(f"Missing RefSeq IDs: {missing_count} ({missing_percentage:.2f}%)")
    
    # Calculate percentage from total RefSeq entries
    found_from_total_percentage = (len(refseq_ids_found) / total_refseq_entries) * 100 if total_refseq_entries > 0 else 0
    print(f"Found RefSeq IDs represent {found_from_total_percentage:.2f}% of total RefSeq entries in the annotation file")
    
    # Print selected TSS count
    print(f"Final TSS position count: {len(tss_positions)}")
    
    # Print group counts and percentages
    print("\nTSS positions found per group:")
    for group in sorted(group_counts.keys(), key=lambda x: (0 if isinstance(x, str) else 1, x if isinstance(x, str) else x)):
        # Percentage of total TSS positions
        group_percentage = (group_counts[group] / len(tss_positions)) * 100 if len(tss_positions) > 0 else 0
        
        # Compare with original count
        original_count = original_group_counts.get(group, 0)
        
        print(f"  Group {group}: {group_counts[group]} TSS positions ({group_percentage:.2f}% of total)")
    
    print(f"\nTotal: {len(tss_positions)} TSS positions")
    
    return tss_positions

def create_tss_bed_files(tss_positions, output_dir, output_bed, original_group_counts, verbose=False):
    """Create BED files with TSS positions and groups"""
    print("\nCreating TSS BED files with expression groups...")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Track unique groups
    unique_groups = set()
    group_counts = defaultdict(int)
    
    # Count TSS positions per group first
    for pos in tss_positions:
        group = pos[6]  # Group is in the 7th column (index 6)
        unique_groups.add(group)
        group_counts[group] += 1
    
    # Debug info about the groups
    if verbose:
        print(f"Groups found in TSS positions: {sorted(unique_groups, key=lambda x: (0 if isinstance(x, str) else 1, x if isinstance(x, str) else x))}")
    
    # Write main BED file with all TSS positions
    output_bed_path = os.path.join(output_dir, output_bed)
    with open(output_bed_path, "w") as f:
        for pos in tss_positions:
            f.write("\t".join(map(str, pos)) + "\n")
    
    print(f"Created main BED file: {output_bed_path}")
    
    # Create separate BED files for each group
    for group in unique_groups:
        group_str = str(group)
        group_file = os.path.join(output_dir, f"group_{group_str}_tss.bed")
        
        # Track unique RefSeq IDs in this BED file but only if verbose is enabled
        if verbose:
            unique_ids = set()
        
        with open(group_file, "w") as f:
            for pos in tss_positions:
                if pos[6] == group:  # Check if this position is in the current group
                    # Include the first 6 columns plus the group name as the last column
                    # Note: This adds group name directly to last column without a header
                    bed_entry = list(map(str, pos[:6]))
                    bed_entry.append(f"Group_{group_str}")
                    f.write("\t".join(bed_entry) + "\n")
                    
                    # Only track unique IDs if verbose is enabled
                    if verbose:
                        refseq_id = pos[3].split('_')[0]
                        unique_ids.add(refseq_id)
        
        # Calculate percentage of total
        percentage = (group_counts[group] / len(tss_positions)) * 100 if len(tss_positions) > 0 else 0
        
        # Print simplified output without detailed duplicate information
        print(f"Created group {group} BED file: {group_file} - {group_counts[group]} TSS positions ({percentage:.2f}% of total)")
    
    # Create a sample file that lists all group BED files
    sample_file = os.path.join(output_dir, "group_samples.txt")
    with open(sample_file, "w") as f:
        # Only write numeric groups in order (skip string groups like "All")
        for group in sorted([g for g in unique_groups if isinstance(g, int)]):
            group_str = str(group)
            f.write(f"{output_dir}/group_{group_str}_tss.bed\t{group_str}\n")
    
    print(f"Created sample file: {sample_file}")
    
    # Print summary with simplified output
    print("\nSummary of TSS BED files by group:")
    total_positions = len(tss_positions)
    for group in sorted(group_counts.keys(), key=lambda x: (0 if isinstance(x, str) else 1, x if isinstance(x, str) else x)):
        percentage = (group_counts[group] / total_positions) * 100 if total_positions > 0 else 0
        print(f"  Group {group}: {group_counts[group]} TSS positions ({percentage:.2f}% of total)")
    
    print(f"\nTotal: {total_positions} TSS positions across {len(unique_groups)} groups")
    
    return output_bed_path

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Print the TSS selection method being used
    print(f"TSS selection method: {args.tss_selection}")
    
    # Get RefSeq gene file - either download or use provided
    refseq_file = args.refseq
    if not refseq_file:
        refseq_file = download_refseq_genes(args.output_dir, args.genome)
    elif not os.path.exists(refseq_file):
        print(f"Error: Provided RefSeq file {refseq_file} does not exist")
        sys.exit(1)
    
    # Read expression group files
    try:
        gene_to_group, original_group_counts = read_expression_group_files(args.expr_dir, args.verbose)
    except Exception as e:
        print(f"Error reading expression group files: {e}")
        sys.exit(1)
    
    # Parse RefSeq genes and extract TSS positions
    try:
        tss_positions = parse_refseq_genes(refseq_file, gene_to_group, original_group_counts, args)
    except Exception as e:
        print(f"Error parsing RefSeq genes: {e}")
        sys.exit(1)
    
    # Create TSS BED files
    try:
        bed_file = create_tss_bed_files(tss_positions, args.output_dir, args.output_bed, original_group_counts, args.verbose)
        if bed_file:
            print(f"\nSuccessfully created RefSeq-based TSS BED files in {args.output_dir}")
            if args.tss_selection == 'all':
                print(f"Note: All alternative TSS positions were retained (using 'all' selection method)")
            else:
                print(f"Note: Used '{args.tss_selection}' method to select representative TSS positions")
        else:
            print("Failed to create TSS BED files.")
            sys.exit(1)
    except Exception as e:
        print(f"Error creating TSS BED files: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
