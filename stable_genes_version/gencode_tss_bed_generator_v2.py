#!/usr/bin/env python3
"""
Script to create TSS BED files from Gencode gene annotations (ENSG IDs) and expression group files.
This script:
1. Downloads Gencode gene annotations if not available
2. Reads expression group files containing ENSG IDs
3. Generates a TSS BED file with expression groups
4. Outputs individual BED files for each expression group
5. Optionally creates a matching RefSeq TSS BED file using gene symbols as mapping

Advanced features include:
- Multiple TSS selection methods for handling alternative transcripts
- Support for different genome builds
- Detailed statistics reporting
- Optional matching RefSeq BED generation
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
    parser = argparse.ArgumentParser(description='Generate TSS BED files from Gencode annotations and expression groups')
    parser.add_argument('--gencode', type=str, help='Path to Gencode gene annotation file (will download if not provided)')
    parser.add_argument('--refseq', type=str, help='Path to RefSeq gene annotation file (for optional RefSeq matching)')
    parser.add_argument('--expr_dir', type=str, required=True, help='Directory containing expression group files')
    parser.add_argument('--output_dir', type=str, default='results', help='Directory to save output BED files')
    parser.add_argument('--output_bed', type=str, default='gencode_tss_with_groups.bed', help='Filename for the main output BED file')
    parser.add_argument('--genome', type=str, default='hg38', help='Genome build (default: hg38)')
    parser.add_argument('--gencode_version', type=str, default='v41', help='Gencode version (default: v41)')
    parser.add_argument('--verbose', action='store_true', help='Print detailed debug information')
    parser.add_argument('--tss_selection', type=str, 
                      choices=['canonical', 'most_upstream', 'most_downstream', 'longest_transcript', 
                               'shortest_transcript', 'first', 'all'], 
                      default='canonical', 
                      help='Method to select representative TSS: canonical (default), most_upstream, most_downstream, '
                           'longest_transcript, shortest_transcript, first, or all (retain all TSS positions)')
    parser.add_argument('--match_refseq', action='store_true', 
                        help='Generate a matching RefSeq BED file for genes found in the Gencode dataset')
    parser.add_argument('--gene_prefix', type=str, default='ENSG',
                        help='Prefix for gene IDs in expression group files (default: ENSG)')
    parser.add_argument('--file_pattern', type=str, default='ENS',
                        help='Pattern in filename to identify group files (default: ENS)')
    return parser.parse_args()

def download_gencode_genes(output_dir, genome="hg38", version="v41"):
    """Download Gencode gene annotation file if not already present, save in a subdirectory"""
    # Create gencode_genefiles subdirectory
    gencode_dir = os.path.join(output_dir, "gencode_genefiles")
    os.makedirs(gencode_dir, exist_ok=True)
    
    # Convert version format if needed (handle both numerical and v-prefixed versions)
    numeric_version = version.lstrip('v')
    v_version = f"v{numeric_version}" if not version.startswith('v') else version
    
    # Determine the correct URL and filename based on genome build
    if genome.lower() in ["hg38", "grch38"]:
        # Standard GRCh38 (hg38) URL
        gencode_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{numeric_version}/gencode.{v_version}.annotation.gtf.gz"
        gencode_file = os.path.join(gencode_dir, f"gencode.{v_version}.{genome}.gtf.gz")
    elif genome.lower() in ["hg19", "grch37"]:
        # GRCh37 (hg19) mapped URL
        gencode_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{numeric_version}/GRCh37_mapping/gencode.{v_version}lift37.annotation.gtf.gz"
        gencode_file = os.path.join(gencode_dir, f"gencode.{v_version}.{genome}.gtf.gz")
    else:
        raise ValueError(f"Unsupported genome build: {genome}")
    
    if not os.path.exists(gencode_file):
        print(f"Downloading Gencode gene annotations for {genome} (version {v_version})...")
        print(f"URL: {gencode_url}")
        download_cmd = f"wget -O {gencode_file} {gencode_url}"
        try:
            subprocess.run(download_cmd, shell=True, check=True)
            print(f"Downloaded Gencode annotations to {gencode_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading Gencode annotations: {e}")
            sys.exit(1)
    else:
        print(f"Using existing Gencode annotations: {gencode_file}")
    
    return gencode_file

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

def read_expression_group_files(expr_dir, file_pattern="ENS", gene_prefix="ENSG", verbose=False):
    """
    Read expression group files and create a mapping from ENSG ID to group.
    Each file should be named with its group number (e.g., "T47D_RNAseq_GeneName_ENS_0.txt").
    Files contain Gencode IDs (ENSG*), one per line.
    """
    print(f"Reading expression group files from {expr_dir}...")
    
    gene_to_group = {}
    group_stats = defaultdict(int)
    groups_found = []
    total_gene_count = 0
    
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
            # Match files like T47D_RNAseq_GeneName_ENS_0.txt
            match = re.search(rf'{file_pattern}_(\d+)\.txt$', filename)
            if match:
                group = int(match.group(1))  # Make sure group is stored as int
        
        if group is None:
            continue
            
        groups_found.append(group)
        count = 0
        
        with open(filepath, 'r') as f:
            for line in f:
                gene_id = line.strip()
                if gene_id:
                    # Ensure gene has correct prefix
                    if not gene_id.startswith(gene_prefix):
                        if verbose:
                            print(f"Warning: Gene ID {gene_id} does not start with expected prefix {gene_prefix}")
                    
                    # Store multiple groups per Gene ID if needed (handle duplicates)
                    gene_to_group[gene_id] = group
                    group_stats[group] += 1
                    count += 1
                    total_gene_count += 1
        
        print(f"Reading group {group} file: {filename} - {count} Gene IDs")
    
    # Count unique Gene IDs
    unique_gene_ids = len(gene_to_group)
    
    print(f"\nFound {len(groups_found)} expression groups: {sorted(groups_found)}")
    print(f"Read {total_gene_count} total Gene IDs from files")
    print(f"Unique Gene IDs: {unique_gene_ids}")
    
    if verbose:
        print("\nGene IDs per group:")
        for group in sorted(group_stats.keys()):
            print(f"  Group {group}: {group_stats[group]} Gene IDs")
    
    # Return both the mapping and the original count statistics for each group
    return gene_to_group, group_stats

def parse_gencode_genes(gencode_file, gene_to_group, original_group_counts, args):
    """
    Parse Gencode gene annotations and extract TSS positions for genes in the expression groups.
    Apply the selected method to handle multiple transcripts for the same ENSG ID.
    Gencode GTF format: chrom, source, feature, start, end, score, strand, frame, attributes
    """
    print(f"Parsing Gencode gene annotations from {gencode_file}...")
    verbose = args.verbose
    
    # List to store TSS positions
    tss_positions = []
    gene_symbols = {}  # Map ENSG IDs to gene symbols
    gene_ids_found = set()
    group_counts = defaultdict(int)
    
    # Store complete entries for selection
    gene_entries = {}
    
    # Open file (handle gzipped files)
    opener = gzip.open if gencode_file.endswith('.gz') else open
    
    # First pass - count transcript entries and genes
    total_gene_count = 0
    total_transcript_count = 0
    
    print("Counting genes and transcripts in annotation file...")
    with opener(gencode_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:  # Ensure we have all necessary fields
                continue
                
            feature_type = fields[2]
            if feature_type == 'gene':
                total_gene_count += 1
            elif feature_type == 'transcript':
                total_transcript_count += 1
    
    print(f"Total genes in annotation file: {total_gene_count}")
    print(f"Total transcripts in annotation file: {total_transcript_count}")
    
    # Second pass - parse and store transcript entries
    print("Parsing gene annotations...")
    with opener(gencode_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            try:
                chrom = fields[0]
                feature_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]
                
                # Skip non-standard chromosomes
                if '_' in chrom or chrom.startswith('chrUn'):
                    continue
                
                # We're only interested in transcript entries
                if feature_type != 'transcript':
                    continue
                
                # Extract gene_id and transcript_id from attributes
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
                transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
                gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
                transcript_type_match = re.search(r'transcript_type "([^"]+)"', attributes)
                
                if not (gene_id_match and transcript_id_match):
                    continue
                    
                gene_id = gene_id_match.group(1)
                transcript_id = transcript_id_match.group(1)
                gene_name = gene_name_match.group(1) if gene_name_match else "Unknown"
                transcript_type = transcript_type_match.group(1) if transcript_type_match else "Unknown"
                
                # Extract version from gene_id if present
                gene_id_no_version = gene_id.split('.')[0]
                
                # Store gene symbol mapping (for later use in RefSeq matching)
                gene_symbols[gene_id_no_version] = gene_name
                
                # Check if this gene ID is in our expression groups
                if gene_id_no_version in gene_to_group:
                    # TSS is at start for + strand genes, end for - strand genes
                    tss = start if strand == "+" else end
                    
                    # Store the group EXACTLY as it appears in the dictionary
                    # This is critical for maintaining type (int vs. str)
                    group = gene_to_group[gene_id_no_version]
                    
                    # Create a key for the gene ID + group combo
                    key = (gene_id_no_version, group)
                    
                    # Calculate transcript length
                    transcript_length = end - start + 1
                    
                    # Create a complete entry for selection
                    entry = {
                        'gene_id': gene_id_no_version,
                        'transcript_id': transcript_id,
                        'gene_name': gene_name,
                        'group': group,
                        'chrom': chrom,
                        'strand': strand,
                        'start': start,
                        'end': end,
                        'tss': tss,
                        'transcript_length': transcript_length,
                        'transcript_type': transcript_type,
                        'is_canonical': transcript_type in ['protein_coding', 'mRNA'],
                        'tss_position': [
                            chrom, tss-1, tss, f"{gene_id_no_version}_{gene_name}", 0, strand, group
                        ]
                    }
                    
                    # Add to entries
                    if key not in gene_entries:
                        gene_entries[key] = []
                    gene_entries[key].append(entry)
                    
                    # Count this as found
                    gene_ids_found.add(gene_id_no_version)
                    
            except Exception as e:
                if verbose:
                    print(f"Error parsing line: {e}")
                continue
    
    # Process entries based on command line args
    if args.tss_selection == 'all':
        print("\nRetaining all TSS positions as requested ('all' selection method)...")
        # Add all TSS positions
        for key, entries in gene_entries.items():
            gene_id, group = key
            for entry in entries:
                tss_positions.append(entry['tss_position'])
                group_counts[group] += 1
            if verbose and len(entries) > 1:
                print(f"  Retained all {len(entries)} TSS positions for {gene_id} (Group {group})")
    else:
        # Select representative TSS positions
        print(f"\nSelecting representative TSS positions using '{args.tss_selection}' method...")
        
        for key, entries in gene_entries.items():
            gene_id, group = key
            
            # If there's only one entry, use it
            if len(entries) == 1:
                tss_positions.append(entries[0]['tss_position'])
                group_counts[group] += 1
                continue
            
            # For multiple entries, apply selection strategy
            selected_entry = None
            
            if args.tss_selection == 'canonical':
                # Prefer protein-coding transcripts
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
                print(f"  Selected 1 of {len(entries)} TSS positions for {gene_id} (Group {group})")
    
    # Calculate percentage of Gene IDs found
    found_percentage = (len(gene_ids_found) / len(gene_to_group)) * 100 if len(gene_to_group) > 0 else 0
    missing_count = len(gene_to_group) - len(gene_ids_found)
    missing_percentage = (missing_count / len(gene_to_group)) * 100 if len(gene_to_group) > 0 else 0
    
    print(f"Found {len(gene_ids_found)} out of {len(gene_to_group)} Gene IDs in the annotations ({found_percentage:.2f}%)")
    print(f"Missing Gene IDs: {missing_count} ({missing_percentage:.2f}%)")
    
    # Print selected TSS count
    print(f"Final TSS position count: {len(tss_positions)}")
    
    # Print group counts and percentages
    print("\nTSS positions found per group:")
    for group in sorted(group_counts.keys()):
        # Percentage of total TSS positions
        group_percentage = (group_counts[group] / len(tss_positions)) * 100 if len(tss_positions) > 0 else 0
        
        # Compare with original count
        original_count = original_group_counts.get(group, 0)
        
        print(f"  Group {group}: {group_counts[group]} TSS positions ({group_percentage:.2f}% of total)")
    
    print(f"\nTotal: {len(tss_positions)} TSS positions")
    
    return tss_positions, gene_symbols, gene_ids_found

def parse_refseq_genes_for_matching(refseq_file, gene_symbols, gene_ids_found, args):
    """
    Parse RefSeq gene annotations and match them with Gencode genes using gene symbols.
    Extract TSS positions for matching genes.
    RefSeq format: bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, etc.
    """
    print(f"Parsing RefSeq gene annotations for matching with Gencode genes...")
    verbose = args.verbose
    
    # Create reverse mapping from gene_symbol to ENSG ID
    symbol_to_ensg = {}
    for ensg_id in gene_ids_found:
        if ensg_id in gene_symbols:
            symbol = gene_symbols[ensg_id]
            symbol_to_ensg[symbol] = ensg_id
    
    print(f"Created mapping for {len(symbol_to_ensg)} gene symbols to ENSG IDs")
    
    # List to store matched TSS positions
    tss_positions = []
    matched_symbols = set()
    
    # Track entries by gene symbol for selection
    symbol_entries = defaultdict(list)
    
    # Open file (handle gzipped files)
    opener = gzip.open if refseq_file.endswith('.gz') else open
    
    print("Extracting RefSeq entries that match Gencode genes...")
    with opener(refseq_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            
            # Check if it's a valid RefSeq line
            if len(fields) < 13:  # Need name, chrom, strand, txStart, txEnd, name2
                continue
                
            try:
                name = fields[1]  # RefSeq ID (NM_* or NR_*)
                chrom = fields[2]
                strand = fields[3]
                tx_start = int(fields[4])
                tx_end = int(fields[5])
                cds_start = int(fields[6])
                cds_end = int(fields[7])
                gene_symbol = fields[12]  # name2 field is the gene symbol
                
                # Skip non-standard chromosomes
                if '_' in chrom or chrom.startswith('chrUn'):
                    continue
                
                # Check if this gene symbol is in our Gencode mapping
                if gene_symbol in symbol_to_ensg:
                    # TSS is at txStart for + strand genes, txEnd for - strand genes
                    tss = tx_start if strand == "+" else tx_end
                    
                    # Get the corresponding ENSG ID
                    ensg_id = symbol_to_ensg[gene_symbol]
                    
                    # Calculate transcript length
                    transcript_length = tx_end - tx_start
                    
                    # Create entry for selection
                    entry = {
                        'refseq_id': name,
                        'gene_symbol': gene_symbol,
                        'ensg_id': ensg_id,
                        'chrom': chrom,
                        'strand': strand,
                        'tx_start': tx_start,
                        'tx_end': tx_end,
                        'tss': tss,
                        'transcript_length': transcript_length,
                        'is_canonical': cds_start < cds_end and cds_start > 0 and cds_end > 0,  # Has valid CDS
                        'tss_position': [
                            chrom, tss-1, tss, f"{name}_{gene_symbol}", 0, strand
                        ]
                    }
                    
                    # Add to entries for this gene symbol
                    symbol_entries[gene_symbol].append(entry)
                    matched_symbols.add(gene_symbol)
                    
            except Exception as e:
                if verbose:
                    print(f"Error parsing RefSeq line: {e}")
                continue
    
    # Select representative TSS for each gene symbol
    print(f"Selecting representative RefSeq TSS positions using '{args.tss_selection}' method...")
    
    for gene_symbol, entries in symbol_entries.items():
        # If there's only one entry, use it
        if len(entries) == 1:
            # Add ENSG ID to the BED entry
            tss_position = entries[0]['tss_position'] + [entries[0]['ensg_id']]
            tss_positions.append(tss_position)
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
        
        # Add the selected entry with ENSG ID
        tss_position = selected_entry['tss_position'] + [selected_entry['ensg_id']]
        tss_positions.append(tss_position)
        
        if verbose and len(entries) > 1:
            print(f"  Selected 1 of {len(entries)} RefSeq TSS positions for {gene_symbol} (ENSG: {selected_entry['ensg_id']})")
    
    print(f"Matched {len(matched_symbols)} gene symbols with RefSeq entries")
    print(f"Created {len(tss_positions)} RefSeq TSS positions for matched Gencode genes")
    
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
        print(f"Groups found in TSS positions: {sorted(unique_groups)}")
    
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
        
        with open(group_file, "w") as f:
            for pos in tss_positions:
                if pos[6] == group:  # Check if this position is in the current group
                    # Include the first 6 columns plus the group name as the last column
                    # Note: This adds group name directly to last column without a header
                    bed_entry = list(map(str, pos[:6]))
                    bed_entry.append(f"Group_{group_str}")
                    f.write("\t".join(bed_entry) + "\n")
        
        # Calculate percentage of total
        percentage = (group_counts[group] / len(tss_positions)) * 100 if len(tss_positions) > 0 else 0
        
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
    for group in sorted(group_counts.keys()):
        percentage = (group_counts[group] / total_positions) * 100 if total_positions > 0 else 0
        print(f"  Group {group}: {group_counts[group]} TSS positions ({percentage:.2f}% of total)")
    
    print(f"\nTotal: {total_positions} TSS positions across {len(unique_groups)} groups")
    
    return output_bed_path

def create_matching_refseq_bed(refseq_tss_positions, tss_positions, output_dir, verbose=False):
    """Create BED file with matching RefSeq TSS positions for Gencode genes"""
    print("\nCreating matching RefSeq BED file...")
    stats_log = []  # Store statistics for log file
    stats_log.append("MATCHING REFSEQ TSS STATISTICS")
    stats_log.append("==============================")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Write matched RefSeq BED file (without groups)
    output_bed_path = os.path.join(output_dir, "refseq_matched_to_gencode.bed")
    with open(output_bed_path, "w") as f:
        for pos in refseq_tss_positions:
            f.write("\t".join(map(str, pos)) + "\n")
    
    msg = f"Created matched RefSeq BED file: {output_bed_path}"
    print(msg)
    stats_log.append(msg)
    
    msg = f"Contains {len(refseq_tss_positions)} TSS positions matched to Gencode genes"
    print(msg)
    stats_log.append(msg)
    
    # Create a mapping file showing the relationship
    mapping_file = os.path.join(output_dir, "refseq_to_gencode_mapping.tsv")
    with open(mapping_file, "w") as f:
        f.write("RefSeq_ID\tGene_Symbol\tENSG_ID\tChromosome\tTSS_Position\tStrand\n")
        for pos in refseq_tss_positions:
            refseq_id = pos[3].split('_')[0]
            gene_symbol = pos[3].split('_', 1)[1] if '_' in pos[3] else "Unknown"
            ensg_id = pos[6]
            chrom = pos[0]
            tss = pos[2]  # The actual TSS position
            strand = pos[5]
            f.write(f"{refseq_id}\t{gene_symbol}\t{ensg_id}\t{chrom}\t{tss}\t{strand}\n")
    
    msg = f"Created RefSeq to Gencode mapping file: {mapping_file}"
    print(msg)
    stats_log.append(msg)
    
    # Build mapping from ENSG ID to group and create group_tss_positions dictionary
    group_tss_positions = defaultdict(list)
    ensg_to_group = {}
    ensg_to_pos = {}
    
    # First build a mapping of ENSG IDs to their groups and positions
    for pos in tss_positions:
        ensg_id = pos[3].split('_')[0]
        ensg_to_group[ensg_id] = pos[6]
        ensg_to_pos[ensg_id] = pos
    
    # Now process RefSeq positions and assign them to groups
    for pos in refseq_tss_positions:
        ensg_id = pos[6]
        if ensg_id in ensg_to_group:
            group = ensg_to_group[ensg_id]
            # Create a new position with the group
            new_pos = pos[:6] + [group]
            group_tss_positions[group].append(new_pos)
        elif verbose:
            print(f"Warning: Could not find group for ENSG ID {ensg_id}, using fallback group 0")
            group_tss_positions[0].append(pos[:6] + [0])
    
    # Create separate BED files for each group
    group_counts = defaultdict(int)
    for group, positions in group_tss_positions.items():
        group_str = str(group)
        group_file = os.path.join(output_dir, f"refseq_matched_group_{group_str}_tss.bed")
        
        with open(group_file, "w") as f:
            for pos in positions:
                # Include the first 6 columns plus the group name as the last column
                bed_entry = list(map(str, pos[:6]))
                bed_entry.append(f"Group_{group_str}")
                f.write("\t".join(bed_entry) + "\n")
        
        group_counts[group] = len(positions)
        msg = f"Created RefSeq matched group {group} BED file: {group_file} - {len(positions)} TSS positions"
        print(msg)
        stats_log.append(msg)
    
    # Create a comprehensive BED file with all groups (matched_refseq_tss_with_groups.bed)
    comprehensive_bed_path = os.path.join(output_dir, "matched_refseq_tss_with_groups.bed")
    with open(comprehensive_bed_path, "w") as f:
        for group, positions in group_tss_positions.items():
            for pos in positions:
                # Write the position with its group (all 7 columns)
                f.write("\t".join(map(str, pos)) + "\n")
    
    msg = f"Created comprehensive RefSeq BED file with groups: {comprehensive_bed_path}"
    print(msg)
    stats_log.append(msg)
    
    # Create side-by-side mapping file with both RefSeq and Gencode data
    side_by_side_file = os.path.join(output_dir, "refseq_gencode_side_by_side.tsv")
    with open(side_by_side_file, "w") as f:
        # Write header
        f.write("RefSeq_ID\tRefSeq_Gene_Symbol\tRefSeq_Chrom\tRefSeq_TSS\tRefSeq_Strand\t")
        f.write("ENSG_ID\tGencode_Gene_Symbol\tGencode_Chrom\tGencode_TSS\tGencode_Strand\tGroup\n")
        
        # Write data
        for pos in refseq_tss_positions:
            refseq_id = pos[3].split('_')[0]
            refseq_symbol = pos[3].split('_', 1)[1] if '_' in pos[3] else "Unknown"
            refseq_chrom = pos[0]
            refseq_tss = pos[2]
            refseq_strand = pos[5]
            ensg_id = pos[6]
            
            # Find corresponding Gencode position
            if ensg_id in ensg_to_pos:
                gencode_pos = ensg_to_pos[ensg_id]
                gencode_symbol = gencode_pos[3].split('_', 1)[1] if '_' in gencode_pos[3] else "Unknown"
                gencode_chrom = gencode_pos[0]
                gencode_tss = gencode_pos[2]
                gencode_strand = gencode_pos[5]
                group = gencode_pos[6]
            else:
                gencode_symbol = "Not_Found"
                gencode_chrom = "Not_Found"
                gencode_tss = "Not_Found"
                gencode_strand = "Not_Found"
                group = "Not_Found"
            
            # Write row
            f.write(f"{refseq_id}\t{refseq_symbol}\t{refseq_chrom}\t{refseq_tss}\t{refseq_strand}\t")
            f.write(f"{ensg_id}\t{gencode_symbol}\t{gencode_chrom}\t{gencode_tss}\t{gencode_strand}\t{group}\n")
    
    msg = f"Created side-by-side RefSeq-Gencode mapping file: {side_by_side_file}"
    print(msg)
    stats_log.append(msg)
    
    # Write statistics to log file
    stats_log_file = os.path.join(output_dir, "refseq_matching_stats.log")
    with open(stats_log_file, "w") as f:
        for line in stats_log:
            f.write(line + "\n")
        
        f.write("\n\nRefSeq TSS positions by group:\n")
        for group in sorted(group_counts.keys()):
            percentage = (group_counts[group] / len(refseq_tss_positions)) * 100 if len(refseq_tss_positions) > 0 else 0
            f.write(f"  Group {group}: {group_counts[group]} positions ({percentage:.2f}% of total)\n")
    
    msg = f"Saved statistics to log file: {stats_log_file}"
    print(msg)
    
    return comprehensive_bed_path

def main():
    """Main function to orchestrate the TSS BED file generation process"""
    # Initialize stats list for logging
    stats_log = []
    
    # Parse command line arguments
    args = parse_args()
    
    # Log start time
    start_time = subprocess.check_output(['date', '+%Y-%m-%d %H:%M:%S']).decode().strip()
    msg = f"Starting TSS BED file generation at {start_time}"
    print(msg)
    stats_log.append(msg)
    stats_log.append(f"TSS selection method: {args.tss_selection}")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Print the TSS selection method being used
    print(f"TSS selection method: {args.tss_selection}")
    
    # Get Gencode gene file - either download or use provided
    gencode_file = args.gencode
    if not gencode_file:
        gencode_file = download_gencode_genes(args.output_dir, args.genome, args.gencode_version)
    elif not os.path.exists(gencode_file):
        print(f"Error: Provided Gencode file {gencode_file} does not exist")
        sys.exit(1)
    
    stats_log.append(f"Using Gencode file: {gencode_file}")
    
    # Read expression group files
    try:
        gene_to_group, original_group_counts = read_expression_group_files(
            args.expr_dir, args.file_pattern, args.gene_prefix, args.verbose
        )
        
        stats_log.append(f"\nRead expression group files from: {args.expr_dir}")
        stats_log.append(f"Found {len(gene_to_group)} unique gene IDs across {len(original_group_counts)} groups")
        
        # Add group stats to log
        stats_log.append("\nGenes per expression group (from input files):")
        for group in sorted(original_group_counts.keys()):
            stats_log.append(f"  Group {group}: {original_group_counts[group]} genes")
        
    except Exception as e:
        print(f"Error reading expression group files: {e}")
        stats_log.append(f"ERROR: Failed to read expression group files: {e}")
        with open(os.path.join(args.output_dir, "gencode_tss_stats.log"), "w") as f:
            for line in stats_log:
                f.write(line + "\n")
        sys.exit(1)
    
    # Parse Gencode genes and extract TSS positions
    try:
        tss_positions, gene_symbols, gene_ids_found = parse_gencode_genes(
            gencode_file, gene_to_group, original_group_counts, args
        )
        
        # Add key statistics to log
        stats_log.append(f"\nFound {len(gene_ids_found)} out of {len(gene_to_group)} Gene IDs in the annotations")
        found_percentage = (len(gene_ids_found) / len(gene_to_group)) * 100 if len(gene_to_group) > 0 else 0
        stats_log.append(f"Gene ID match rate: {found_percentage:.2f}%")
        stats_log.append(f"Final TSS position count: {len(tss_positions)}")
        
        # Calculate group stats for log
        group_counts = defaultdict(int)
        for pos in tss_positions:
            group_counts[pos[6]] += 1
            
        stats_log.append("\nTSS positions found per group:")
        for group in sorted(group_counts.keys()):
            percentage = (group_counts[group] / len(tss_positions)) * 100 if len(tss_positions) > 0 else 0
            stats_log.append(f"  Group {group}: {group_counts[group]} TSS positions ({percentage:.2f}% of total)")
        
    except Exception as e:
        print(f"Error parsing Gencode genes: {e}")
        stats_log.append(f"ERROR: Failed to parse Gencode genes: {e}")
        with open(os.path.join(args.output_dir, "gencode_tss_stats.log"), "w") as f:
            for line in stats_log:
                f.write(line + "\n")
        sys.exit(1)
    
    # Create TSS BED files
    try:
        bed_file = create_tss_bed_files(
            tss_positions, args.output_dir, args.output_bed, original_group_counts, args.verbose
        )
        if bed_file:
            msg = f"\nSuccessfully created Gencode-based TSS BED files in {args.output_dir}"
            print(msg)
            stats_log.append(msg)
            
            if args.tss_selection == 'all':
                msg = f"Note: All alternative TSS positions were retained (using 'all' selection method)"
                print(msg)
                stats_log.append(msg)
            else:
                msg = f"Note: Used '{args.tss_selection}' method to select representative TSS positions"
                print(msg)
                stats_log.append(msg)
                
            stats_log.append(f"Main BED file created: {bed_file}")
        else:
            print("Failed to create TSS BED files.")
            stats_log.append("ERROR: Failed to create TSS BED files")
            with open(os.path.join(args.output_dir, "gencode_tss_stats.log"), "w") as f:
                for line in stats_log:
                    f.write(line + "\n")
            sys.exit(1)
    except Exception as e:
        print(f"Error creating TSS BED files: {e}")
        stats_log.append(f"ERROR: Failed to create TSS BED files: {e}")
        with open(os.path.join(args.output_dir, "gencode_tss_stats.log"), "w") as f:
            for line in stats_log:
                f.write(line + "\n")
        sys.exit(1)
    
    # Optionally create matching RefSeq BED file
    if args.match_refseq:
        refseq_stats_log = []
        refseq_stats_log.append("REFSEQ MATCHING STATISTICS")
        refseq_stats_log.append("=========================")
        
        msg = "\nGenerating matching RefSeq TSS BED file..."
        print(msg)
        stats_log.append(msg)
        refseq_stats_log.append(msg)
        
        # Get RefSeq gene file - either download or use provided
        refseq_file = args.refseq
        if not refseq_file:
            refseq_file = download_refseq_genes(args.output_dir, args.genome)
        elif not os.path.exists(refseq_file):
            print(f"Error: Provided RefSeq file {refseq_file} does not exist")
            stats_log.append(f"ERROR: Provided RefSeq file {refseq_file} does not exist")
            with open(os.path.join(args.output_dir, "gencode_tss_stats.log"), "w") as f:
                for line in stats_log:
                    f.write(line + "\n")
            sys.exit(1)
        
        refseq_stats_log.append(f"Using RefSeq file: {refseq_file}")
        
        try:
            # Parse RefSeq genes and match with Gencode genes
            refseq_tss_positions = parse_refseq_genes_for_matching(
                refseq_file, gene_symbols, gene_ids_found, args
            )
            
            refseq_stats_log.append(f"Found {len(refseq_tss_positions)} RefSeq TSS positions matching Gencode genes")
            match_percentage = (len(refseq_tss_positions) / len(gene_ids_found)) * 100 if len(gene_ids_found) > 0 else 0
            refseq_stats_log.append(f"RefSeq match rate: {match_percentage:.2f}% of found Gencode genes")

            # Get detailed statistics about RefSeq matches
            group_counts = defaultdict(int)
            for pos in refseq_tss_positions:
                ensg_id = pos[6]
                
                # Find corresponding group from the Gencode data
                for orig_pos in tss_positions:
                    if orig_pos[3].split('_')[0] == ensg_id:
                        group = orig_pos[6]
                        group_counts[group] += 1
                        break

            # Add group statistics to the log
            total_positions = len(refseq_tss_positions)
            refseq_stats_log.append("\nRefSeq TSS positions found per group:")
            for group in sorted(group_counts.keys()):
                percentage = (group_counts[group] / total_positions) * 100 if total_positions > 0 else 0
                refseq_stats_log.append(f"  Group {group}: {group_counts[group]} TSS positions ({percentage:.2f}% of total)")


            # Write initial RefSeq stats to log file
            refseq_stats_file = os.path.join(args.output_dir, "refseq_tss_stats.log")
            with open(refseq_stats_file, "w") as f:
                for line in refseq_stats_log:
                    f.write(line + "\n")
                    
            print(f"Saved initial RefSeq matching statistics to: {refseq_stats_file}")
            
            # Create matching RefSeq BED file
            matched_bed = create_matching_refseq_bed(
                refseq_tss_positions, tss_positions, args.output_dir, args.verbose
            )
            
            msg = f"\nSuccessfully created matching RefSeq TSS BED files in {args.output_dir}"
            print(msg)
            stats_log.append(msg)
            refseq_stats_log.append(msg)
            refseq_stats_log.append(f"Main RefSeq matched BED file created: {matched_bed}")
                
            # Write RefSeq stats to log file
            refseq_stats_file = os.path.join(args.output_dir, "refseq_matching_stats.log")
            with open(refseq_stats_file, "w") as f:
                for line in refseq_stats_log:
                    f.write(line + "\n")
                    
            print(f"Saved RefSeq matching statistics to: {refseq_stats_file}")
            
        except Exception as e:
            print(f"Error creating matching RefSeq BED files: {e}")
            stats_log.append(f"ERROR: Failed to create matching RefSeq BED files: {e}")
            with open(os.path.join(args.output_dir, "gencode_tss_stats.log"), "w") as f:
                for line in stats_log:
                    f.write(line + "\n")
            sys.exit(1)
    
    # Write stats log file
    stats_log_file = os.path.join(args.output_dir, "gencode_tss_stats.log")
    with open(stats_log_file, "w") as f:
        for line in stats_log:
            f.write(line + "\n")
    
    print(f"Saved statistics to log file: {stats_log_file}")
    
    # Log end time
    end_time = subprocess.check_output(['date', '+%Y-%m-%d %H:%M:%S']).decode().strip()
    print(f"\nTSS BED file generation completed successfully at {end_time}")


if __name__ == "__main__":
    start_time = subprocess.check_output(['date', '+%Y-%m-%d %H:%M:%S']).decode().strip()
    print(f"Starting TSS BED file generation at {start_time}")
    
    try:
        main()
        
        end_time = subprocess.check_output(['date', '+%Y-%m-%d %H:%M:%S']).decode().strip()
        print(f"\nTSS BED file generation completed successfully at {end_time}")
        print("All files generated successfully!")
    except Exception as e:
        print(f"\nERROR: TSS BED file generation failed: {str(e)}")
        sys.exit(1)

    
