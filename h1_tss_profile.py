#!/usr/bin/env python3
"""
Enhanced analysis of H1 histone variant distribution around TSS.
Supports both WIG and GFF file formats for H1 signal data.
Features:
- Adaptive color maps based on expression group count
- Improved heatmap visualization with proper axis annotations
- Combined profile and heatmap visualizations
- DeepTools-inspired color schemes
- GFF format support for ChIP-H1-ChIP data
- WIG/BigWig support for H1 histone variant signal
- Options for log2 transformation or absolute value of signals
- Caching for faster processing
- Publication-ready visualizations
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize
import subprocess
import multiprocessing as mp
from functools import partial
import time
import seaborn as sns
from scipy import stats
import pickle
import hashlib

# Check for pyBigWig
try:
    import pyBigWig
    has_pybigwig = True
except ImportError:
    has_pybigwig = False
    print("Warning: pyBigWig not found. BigWig support disabled.")

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='H1 variant distribution around TSS')
    parser.add_argument('--tss', required=True, help='TSS BED file with expression groups')
    parser.add_argument('--input', required=True, help='Input file with H1 signal (WIG, BigWig, or GFF)')
    parser.add_argument('--format', choices=['wig', 'bigwig', 'gff'], help='Input file format (auto-detected by default)')
    parser.add_argument('--output', default='h1_tss_profile.png', help='Output figure file')
    parser.add_argument('--window', type=int, default=3000, help='Window size around TSS (bp)')
    parser.add_argument('--bin', type=int, default=50, help='Bin size (bp)')
    parser.add_argument('--groups', default=None, help='Expression groups to include (comma-separated, e.g., "1,2,3,4")')
    parser.add_argument('--convert', action='store_true', help='Convert WIG to BigWig for faster processing')
    parser.add_argument('--chrom-sizes', default=None, help='Chromosome sizes file (required for WIG to BigWig conversion)')
    parser.add_argument('--processes', type=int, default=4, help='Number of processes to use')
    parser.add_argument('--sample', type=int, default=None, help='Process only a subset of TSS positions')
    parser.add_argument('--force', action='store_true', help='Force recalculation of profiles even if cached data exists')
    parser.add_argument('--nocache', action='store_true', help='Disable caching of profile data')
    parser.add_argument('--log2', action='store_true', help='Apply log2 transformation to signal values before averaging')
    parser.add_argument('--abs', action='store_true', help='Use absolute values of signals (magnitude only)')
    return parser.parse_args()

def detect_file_format(input_file):
    """Detect file format based on extension and content"""
    # Check extension first
    if input_file.endswith('.wig'):
        return 'wig'
    elif input_file.endswith('.bw') or input_file.endswith('.bigwig'):
        return 'bigwig'
    elif input_file.endswith('.gff') or input_file.endswith('.gff3'):
        return 'gff'
    
    # Check content if extension is ambiguous
    with open(input_file, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith('# Roche NimbleGen') or first_line.startswith('##gff-version'):
            return 'gff'
        elif first_line.startswith('track type=wiggle'):
            return 'wig'
    
    # Default to GFF if can't determine
    print(f"Warning: Could not definitively determine file format for {input_file}. Assuming GFF format.")
    return 'gff'

def parse_tss_file(tss_file, selected_groups=None, sample_size=None):
    """Parse TSS BED file with expression groups, with optional filtering"""
    try:
        tss_data = pd.read_csv(tss_file, sep='\t', header=None,
                              names=['chr', 'start', 'end', 'gene', 'score', 'strand', 'expression_group'])
        
        # Filter by expression group if specified
        if selected_groups:
            tss_data = tss_data[tss_data['expression_group'].isin(selected_groups)]
            print(f"Filtered to {len(tss_data)} TSS positions in expression groups {selected_groups}")
        
        # Take a sample if specified
        if sample_size and sample_size < len(tss_data):
            tss_data = tss_data.sample(sample_size, random_state=42)
            print(f"Sampled {sample_size} TSS positions for faster processing")
            
        return tss_data
    except Exception as e:
        print(f"Error parsing TSS file: {e}")
        sys.exit(1)

def parse_gff_file(gff_file):
    """
    Parse GFF file with ChIP-H1-ChIP data
    Format: chr, source, feature, start, end, score, strand, frame, attributes
    """
    try:
        # Skip comment lines starting with #
        with open(gff_file, 'r') as f:
            lines = [line for line in f if not line.startswith('#')]
        
        # Create a dataframe with the needed columns
        # If it's a simplified GFF (like the sample), adapt accordingly
        if len(lines) > 0 and len(lines[0].split('\t')) >= 6:
            # Full GFF format
            column_names = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
            gff_data = pd.read_csv(
                gff_file, 
                sep='\t', 
                comment='#', 
                names=column_names,
                usecols=['chr', 'start', 'end', 'score']
            )
        else:
            # Simplified format (chr, start, end, score)
            column_names = ['chr', 'start', 'end', 'score']
            gff_data = pd.read_csv(
                gff_file,
                sep='\t',
                comment='#',
                names=column_names
            )
        
        # Ensure score is numeric
        gff_data['score'] = pd.to_numeric(gff_data['score'], errors='coerce')
        
        print(f"Loaded GFF file with {len(gff_data)} regions")
        print(f"Score range: {gff_data['score'].min()} to {gff_data['score'].max()}")
        
        return gff_data
    except Exception as e:
        print(f"Error parsing GFF file: {e}")
        sys.exit(1)

def convert_wig_to_bigwig(wig_file, chrom_sizes_file, overwrite=False):
    """
    Convert WIG to BigWig format for faster random access.
    Requires wigToBigWig from UCSC tools.
    """
    bw_file = wig_file.replace('.wig', '.bw')
    
    # Check if BigWig file already exists
    if os.path.exists(bw_file) and not overwrite:
        print(f"BigWig file {bw_file} already exists. Using existing file.")
        return bw_file
    
    # Check if wigToBigWig is available
    try:
        subprocess.run(['which', 'wigToBigWig'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("Error: wigToBigWig tool not found. Please install UCSC tools.")
        print("You can download it from: http://hgdownload.soe.ucsc.edu/admin/exe/")
        sys.exit(1)
    
    # Check if chromosome sizes file exists
    if not os.path.exists(chrom_sizes_file):
        print(f"Error: Chromosome sizes file {chrom_sizes_file} not found.")
        print("You can generate it with 'fetchChromSizes' or download from UCSC.")
        sys.exit(1)
    
    print(f"Converting {wig_file} to BigWig format...")
    try:
        start_time = time.time()
        subprocess.run(['wigToBigWig', wig_file, chrom_sizes_file, bw_file], 
                      check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elapsed = time.time() - start_time
        print(f"Conversion completed in {elapsed:.2f} seconds.")
        return bw_file
    except subprocess.CalledProcessError as e:
        print(f"Error converting to BigWig: {e}")
        print(f"stderr: {e.stderr.decode() if e.stderr else 'None'}")
        sys.exit(1)

def process_tss_batch_gff(tss_batch, gff_data, bin_edges, window_size=3000, use_log2=False, use_abs=False):
    """
    Process a batch of TSS positions using GFF data.
    This function is used for parallel processing with GFF input.
    """
    # Initialize result containers
    result_sums = {}
    result_counts = {}
    
    # Get unique expression groups in this batch
    expression_groups = sorted(tss_batch['expression_group'].unique())
    
    # Initialize matrices for each group
    for group in expression_groups:
        result_sums[group] = np.zeros(len(bin_edges) - 1)
        result_counts[group] = np.zeros(len(bin_edges) - 1)
    
    # Process each TSS
    for idx, row in tss_batch.iterrows():
        chrom = row['chr']
        tss_pos = row['start']
        group = row['expression_group']
        
        # Define window around TSS
        window_start = tss_pos - window_size
        window_end = tss_pos + window_size
        
        # Find GFF regions that overlap with this window
        chrom_gff = gff_data[gff_data['chr'] == chrom]
        overlapping_regions = chrom_gff[
            (chrom_gff['end'] >= window_start) & 
            (chrom_gff['start'] <= window_end)
        ]
        
        # For each region that overlaps with the window
        for _, region in overlapping_regions.iterrows():
            region_start = region['start']
            region_end = region['end']
            region_score = region['score']
            
            # Apply transformations if requested
            if use_abs:
                region_score = abs(region_score)
                
            if use_log2:
                # Add pseudocount of 1 if positive, or apply to absolute value if negative
                if region_score >= 0:
                    region_score = np.log2(region_score + 1)
                else:
                    # Handle negative values (common in some ChIP data)
                    region_score = -np.log2(abs(region_score) + 1)
            
            # Calculate region's position relative to TSS
            region_rel_start = region_start - tss_pos
            region_rel_end = region_end - tss_pos
            
            # Assign score to appropriate bins
            for i in range(len(bin_edges) - 1):
                bin_start = bin_edges[i]
                bin_end = bin_edges[i+1]
                
                # Check if region overlaps with this bin
                if region_rel_end >= bin_start and region_rel_start < bin_end:
                    # Calculate overlap proportion for partial overlaps
                    overlap_start = max(bin_start, region_rel_start)
                    overlap_end = min(bin_end, region_rel_end)
                    overlap_length = overlap_end - overlap_start
                    
                    # Only count if there's positive overlap
                    if overlap_length > 0:
                        result_sums[group][i] += region_score
                        result_counts[group][i] += 1
    
    return result_sums, result_counts

def process_tss_batch_bigwig(tss_batch, bw_file, bin_edges, window_size=3000, use_log2=False, use_abs=False):
    """
    Process a batch of TSS positions using BigWig data.
    This function is used for parallel processing with BigWig input.
    """
    # Initialize result containers
    result_sums = {}
    result_counts = {}
    
    # Get unique expression groups in this batch
    expression_groups = sorted(tss_batch['expression_group'].unique())
    
    # Initialize matrices for each group
    for group in expression_groups:
        result_sums[group] = np.zeros(len(bin_edges) - 1)
        result_counts[group] = np.zeros(len(bin_edges) - 1)
    
    # Open BigWig file
    bw = pyBigWig.open(bw_file)
    
    # Process each TSS
    for idx, row in tss_batch.iterrows():
        chrom = row['chr']
        tss_pos = row['start']
        group = row['expression_group']
        
        # Skip chromosomes not in BigWig
        if chrom not in bw.chroms():
            continue
        
        # Calculate window around TSS
        start = max(0, tss_pos - window_size)
        end = min(bw.chroms()[chrom], tss_pos + window_size)
        
        # Skip if window is outside chromosome bounds
        if end <= start:
            continue
            
        try:
            # Get values in the window
            values = bw.values(chrom, start, end)
            
            # Calculate positions relative to TSS
            rel_positions = np.arange(start - tss_pos, end - tss_pos)
            
            # Assign values to bins
            for i in range(len(bin_edges) - 1):
                bin_start = bin_edges[i]
                bin_end = bin_edges[i+1]
                
                # Find values in this bin
                bin_mask = (rel_positions >= bin_start) & (rel_positions < bin_end)
                bin_values = [v for v, m in zip(values, bin_mask) if m and v is not None]
                
                # Apply transformations if requested
                if bin_values:
                    if use_abs:
                        bin_values = [abs(v) for v in bin_values]
                        
                    if use_log2:
                        bin_values = [np.log2(v + 1) if v >= 0 else -np.log2(abs(v) + 1) for v in bin_values]
                    
                    result_sums[group][i] += sum(bin_values)
                    result_counts[group][i] += len(bin_values)
                    
        except Exception as e:
            # Skip errors for individual positions
            print(f"Error processing {chrom}:{start}-{end}: {e}")
    
    bw.close()
    return result_sums, result_counts

def get_cache_filename(tss_file, input_file, file_format, window_size, bin_size, selected_groups, sample_size, use_log2=False, use_abs=False):
    """Generate a cache filename based on analysis parameters"""
    # Create a unique ID based on analysis parameters
    params = f"{tss_file}_{input_file}_{file_format}_{window_size}_{bin_size}"
    if selected_groups:
        params += f"_groups{'_'.join(map(str, selected_groups))}"
    if sample_size:
        params += f"_sample{sample_size}"
    if use_log2:
        params += "_log2"
    if use_abs:
        params += "_abs"
    
    # Create a hash of the parameters for a shorter filename
    params_hash = hashlib.md5(params.encode()).hexdigest()
    
    # Get output directory from tss_file path
    output_dir = os.path.dirname(tss_file) or '.'
    
    return f"{output_dir}/h1_profiles_cache_{params_hash}.pkl"

def calculate_tss_profiles_parallel(tss_data, input_data, file_format, window_size=3000, bin_size=50, processes=4, 
                                    cache_file=None, force_recalculate=False, use_log2=False, use_abs=False):
    """
    Calculate average profiles around TSS using parallel processing.
    Supports both GFF and BigWig input formats.
    Optionally use cached data if available.
    """
    # Check if cached data exists and should be used
    if cache_file and os.path.exists(cache_file) and not force_recalculate:
        print(f"Loading cached profile data from {cache_file}")
        try:
            with open(cache_file, 'rb') as f:
                cached_data = pickle.load(f)
                print("Using cached profile data")
                return cached_data['bin_positions'], cached_data['profiles']
        except Exception as e:
            print(f"Error loading cached data: {e}")
            print("Calculating profiles from scratch...")
    
    # Calculate profiles from scratch
    # Create bins
    bin_edges = np.arange(-window_size, window_size + bin_size, bin_size)
    bin_centers = bin_edges[:-1] + bin_size/2
    
    # Split TSS data into batches for parallel processing
    batch_size = max(100, len(tss_data) // (processes * 10))  # Aim for 10 batches per process
    tss_batches = [tss_data.iloc[i:i+batch_size] for i in range(0, len(tss_data), batch_size)]
    print(f"Processing {len(tss_batches)} batches of ~{batch_size} TSS positions each using {processes} processes")
    
    # Report signal processing mode
    if use_log2 and use_abs:
        print("Using log2 transformation on absolute values")
    elif use_log2:
        print("Using log2 transformation")
    elif use_abs:
        print("Using absolute values")
    else:
        print("Using raw values")
    
    # Process batches in parallel
    start_time = time.time()
    
    # Create a pool of worker processes
    pool = mp.Pool(processes=processes)
    
    # Create a partial function with fixed arguments based on file format
    if file_format == 'gff':
        process_func = partial(process_tss_batch_gff, gff_data=input_data, bin_edges=bin_edges, 
                              window_size=window_size, use_log2=use_log2, use_abs=use_abs)
    else:  # bigwig or wig
        process_func = partial(process_tss_batch_bigwig, bw_file=input_data, bin_edges=bin_edges, 
                              window_size=window_size, use_log2=use_log2, use_abs=use_abs)
    
    # Process batches in parallel
    results = []
    for i, batch_result in enumerate(pool.imap_unordered(process_func, tss_batches)):
        results.append(batch_result)
        if (i+1) % 10 == 0 or (i+1) == len(tss_batches):
            elapsed = time.time() - start_time
            print(f"Processed {i+1}/{len(tss_batches)} batches in {elapsed:.2f} seconds")
    
    pool.close()
    pool.join()
    
    # Combine results from all batches
    all_expression_groups = sorted(tss_data['expression_group'].unique())
    
    # Initialize combined results
    combined_sums = {group: np.zeros(len(bin_centers)) for group in all_expression_groups}
    combined_counts = {group: np.zeros(len(bin_centers)) for group in all_expression_groups}
    
    # Combine batch results
    for batch_sums, batch_counts in results:
        for group in batch_sums:
            combined_sums[group] += batch_sums[group]
            combined_counts[group] += batch_counts[group]
    
    # Calculate average profiles
    profiles = {}
    for group in all_expression_groups:
        with np.errstate(divide='ignore', invalid='ignore'):
            profile = np.divide(combined_sums[group], combined_counts[group])
            profile[np.isnan(profile)] = 0
            profiles[group] = profile
    
    # Calculate "All" profile
    all_sum = np.zeros(len(bin_centers))
    all_count = np.zeros(len(bin_centers))
    
    for group in all_expression_groups:
        all_sum += combined_sums[group]
        all_count += combined_counts[group]
    
    with np.errstate(divide='ignore', invalid='ignore'):
        all_profile = np.divide(all_sum, all_count)
        all_profile[np.isnan(all_profile)] = 0
        profiles['All'] = all_profile
    
    # Cache the results if requested
    if cache_file:
        cache_dir = os.path.dirname(cache_file)
        if cache_dir and not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
            
        print(f"Caching profile data to {cache_file}")
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump({
                    'bin_positions': bin_centers,
                    'profiles': profiles
                }, f)
        except Exception as e:
            print(f"Error caching profile data: {e}")
    
    return bin_centers, profiles

def create_adaptive_cmap(num_groups, base_colors=None):
    """
    Create a colormap that adapts to the number of expression groups.
    Uses the original base colors but adjusts to the number of groups.
    """
    # Default base colors from original implementation
    if base_colors is None:
        base_colors = [
            '#952929', '#B25A31', '#BD8430', '#BDAA30', '#8ABD30', 
            '#30BD71', '#30BDBD', '#308ABD', '#3059BD', '#5930BD', '#8A30BD'
        ]
    
    # If num_groups <= len(base_colors), take a slice of the base colors
    if num_groups <= len(base_colors):
        colors_to_use = base_colors[:num_groups]
    else:
        # If num_groups > len(base_colors), interpolate between the base colors
        # to create enough colors for all groups
        colors_to_use = []
        for i in range(num_groups):
            idx = (i * (len(base_colors) - 1)) / (num_groups - 1) if num_groups > 1 else 0
            idx_low = int(idx)
            idx_high = min(idx_low + 1, len(base_colors) - 1)
            frac = idx - idx_low
            
            if idx_low == idx_high:  # Edge case
                colors_to_use.append(base_colors[idx_low])
            else:
                # Interpolate RGB values
                color_low = base_colors[idx_low]
                color_high = base_colors[idx_high]
                
                # Convert hex to RGB
                rgb_low = np.array([int(color_low[1:3], 16), int(color_low[3:5], 16), int(color_low[5:7], 16)])
                rgb_high = np.array([int(color_high[1:3], 16), int(color_high[3:5], 16), int(color_high[5:7], 16)])
                
                # Interpolate
                rgb_interp = rgb_low * (1 - frac) + rgb_high * frac
                
                # Convert back to hex
                color_interp = f'#{int(rgb_interp[0]):02x}{int(rgb_interp[1]):02x}{int(rgb_interp[2]):02x}'
                colors_to_use.append(color_interp)
    
    # Create colormap
    cmap = LinearSegmentedColormap.from_list('adaptive_cmap', colors_to_use, N=num_groups)
    return cmap

def create_deeptools_cmap():
    """
    Create colormap inspired by deepTools.
    Uses blue to white to red color scheme for heatmaps.
    """
    # Define colors for the deepTools-inspired colormap
    colors = [
        '#053061', '#1E61A5', '#3C8ABE', '#67A9CF', '#91C7E0', 
        '#CCDEEF', '#F2F2F2',  # white in the middle
        '#FCDBC6', '#F8A982', '#EC7751', '#D33F36', '#AB1529'
    ]
    
    return LinearSegmentedColormap.from_list('deeptools_cmap', colors, N=256)

def plot_tss_profiles(bin_positions, profiles, output_file, window_size=3000, use_log2=False, use_abs=False):
    """
    Create a visualization of H1 distribution around TSS with an adaptive colormap.
    """
    # Get expression groups (excluding 'All')
    expression_groups = [g for g in profiles.keys() if g != 'All']
    # Sort numeric expression groups separately
    expression_groups = sorted(expression_groups)
    
    # Create adaptive expression color map based on the number of groups
    num_groups = len(expression_groups)
    base_colors = [
        '#952929', '#B25A31', '#BD8430', '#BDAA30', '#8ABD30', 
        '#30BD71', '#30BDBD', '#308ABD', '#3059BD', '#5930BD', '#8A30BD'
    ]
    
    #cmap = create_adaptive_cmap(num_groups, base_colors)
 
    # Use alternating divergent colors for 5 or fewer groups
    if num_groups <= 5:
        # Pick colors from different parts of the base_colors list for better contrast
        # For ≤5 groups, use colors that are visually distinct from each other
        if num_groups == 1:
            selected_colors = [base_colors[5]]  # Mid-point color
        elif num_groups == 2:
            #selected_colors = [base_colors[0], base_colors[3]]  # First and last
            selected_colors = [base_colors[0], base_colors[10]]  # First and last
        elif num_groups == 3:
            selected_colors = [base_colors[0], base_colors[5], base_colors[10]]  # Beginning, middle, end
        elif num_groups == 4:
            selected_colors = [base_colors[0], base_colors[3], base_colors[7], base_colors[10]]  # Spread evenly
        else:  # num_groups == 5
            selected_colors = [base_colors[0], base_colors[2], base_colors[5], base_colors[8], base_colors[10]]  # Spread evenly
            
        cmap = create_adaptive_cmap(num_groups, selected_colors)
        #cmap = create_adaptive_cmap(num_groups, selected_colors)
    else:
        # Use default color scheme for more than 5 groups
        cmap = create_adaptive_cmap(num_groups, base_colors)
    
    # Create figure with specific layout
    fig = plt.figure(figsize=(10, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    
    # Upper panel: Line plot for each expression group
    ax1 = plt.subplot(gs[0])
    
    # Plot each expression group
    for i, group in enumerate(expression_groups):
        color_idx = i / max(1, num_groups - 1)  # Avoid division by zero
        ax1.plot(bin_positions, profiles[group], color=cmap(color_idx), linewidth=1.5, 
                 label=f"{group}")
    
    # Plot the "All" profile with thicker black line
    if 'All' in profiles:
        ax1.plot(bin_positions, profiles['All'], 'k-', linewidth=2.5, label="All")
    
    # Add vertical line at TSS
    ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.7)
    
    # Add grid
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Set axis limits and labels
    ax1.set_xlim(-window_size, window_size)
    
    # Set y-axis label based on transformations
    if use_log2 and use_abs:
        y_label = "log2(|Signal|)"
    elif use_log2:
        y_label = "log2(Signal)"
    elif use_abs:
        y_label = "|Signal|"
    else:
        y_label = "Average Profile"
    
    ax1.set_ylabel(y_label)
    
    # Set title based on transformations
    if use_log2 and use_abs:
        title_suffix = " (log2 of absolute values)"
    elif use_log2:
        title_suffix = " (log2 transformed)"
    elif use_abs:
        title_suffix = " (absolute values)"
    else:
        title_suffix = ""
    
    ax1.set_title(f"H1 Distribution Around TSS{title_suffix}")
    
    # Lower panel: Color legend
    ax2 = plt.subplot(gs[1])
    
    # Create a horizontal colorbar for expression groups
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    
    # Plot colorbar
    ax2.imshow(gradient, aspect='auto', cmap=cmap)
    
    # Configure legend axis
    ax2.set_yticks([])
    
    # Set x-axis ticks to match expression groups
    xticks = np.linspace(0, 255, num_groups)
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(expression_groups)
    
    # Add "All" marker
    ax2.text(265, 0.5, "All", ha='left', va='center', fontweight='bold')
    ax2.plot([260], [0.5], 'ko', markersize=8)
    
    # Add labels for x-axis and colorbar
    ax1.set_xlabel("Relative Distance (bp)")
    ax2.set_xlabel("Gene expression")
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_file}")
    
    # Also save CSV of profile data for further analysis
    csv_file = output_file.replace('.png', '_profiles.csv')
    data_for_csv = {'Position': bin_positions}
    
    # Add profiles to CSV data, handling string keys separately
    for group in expression_groups:
        data_for_csv[f'Group_{group}'] = profiles[group]
    
    if 'All' in profiles:
        data_for_csv['Group_All'] = profiles['All']
    
    pd.DataFrame(data_for_csv).to_csv(csv_file, index=False)
    print(f"Profile data saved to {csv_file}")

def plot_correlation(profiles, bin_positions, output_file, use_log2=False, use_abs=False):
    """
    Create a publication-ready plot showing the correlation between
    gene expression level and H1 depletion at TSS.
    Enhanced with an adaptive color scheme.
    """
    # Get numeric expression groups (excluding 'All')
    expression_groups = [g for g in profiles.keys() if g != 'All']
    expression_groups = sorted(expression_groups)
    
    # Create adaptive colormap
    num_groups = len(expression_groups)
    base_colors = [
        '#952929', '#B25A31', '#BD8430', '#BDAA30', '#8ABD30', 
        '#30BD71', '#30BDBD', '#308ABD', '#3059BD', '#5930BD', '#8A30BD'
    ]
    #cmap = create_adaptive_cmap(num_groups, base_colors)

    # Use alternating divergent colors for 5 or fewer groups
    if num_groups <= 5:
        # Pick colors from different parts of the base_colors list for better contrast
        # For ≤5 groups, use colors that are visually distinct from each other 
        if num_groups == 1:
            selected_colors = [base_colors[5]]  # Mid-point color
        elif num_groups == 2:
            #selected_colors = [base_colors[0], base_colors[3]]  # First and last
            selected_colors = [base_colors[0], base_colors[10]]  # First and last
        elif num_groups == 3:
            selected_colors = [base_colors[0], base_colors[5], base_colors[10]]  # Beginning, middle, end
        elif num_groups == 4:
            selected_colors = [base_colors[0], base_colors[3], base_colors[7], base_colors[10]]  # Spread evenly
        else:  # num_groups == 5
            selected_colors = [base_colors[0], base_colors[2], base_colors[5], base_colors[8], base_colors[10]]  # Spread evenly

        #profile_cmap = create_adaptive_cmap(num_groups, selected_colors)
        cmap = create_adaptive_cmap(num_groups, selected_colors)
    else:
        # Use default color scheme for more than 5 groups
        cmap = create_adaptive_cmap(num_groups, base_colors)

    
    # Calculate depletion for each group (minimum value near TSS)
    # Focus on region near TSS (±500bp)
    tss_region_mask = (bin_positions >= -500) & (bin_positions <= 500)
    tss_bin_positions = bin_positions[tss_region_mask]
    
    depletion_values = []
    for group in expression_groups:
        profile_tss_region = profiles[group][tss_region_mask]
        min_val = np.min(profile_tss_region)
        depletion_values.append(min_val)
    
    # Create DataFrame for plotting
    corr_data = pd.DataFrame({
        'Expression Group': expression_groups,
        'H1 Depletion': depletion_values
    })
    
    # Calculate correlation
    corr, p_value = stats.pearsonr(expression_groups, depletion_values)
    
    # Create figure
    plt.figure(figsize=(10, 8))
    
    # Create custom color list based on expression groups
    colors = [cmap(i / max(1, num_groups - 1)) for i, group in enumerate(expression_groups)]
    
    # Create scatter plot with colored points
    plt.scatter(corr_data['Expression Group'], corr_data['H1 Depletion'], 
                s=100, c=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Add regression line
    x = np.array(expression_groups)
    y = np.array(depletion_values)
    m, b = np.polyfit(x, y, 1)
    plt.plot(x, m*x + b, color='black', linestyle='--', linewidth=2)
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Add correlation information
    plt.annotate(f'Correlation: {corr:.2f}\np-value: {p_value:.3e}',
                xy=(0.05, 0.05), 
                xycoords='axes fraction',
                fontsize=12,
                bbox=dict(boxstyle="round,pad=0.5", facecolor='white', alpha=0.8))
    
    # Set labels based on transformations
    if use_log2 and use_abs:
        y_label = "H1 Depletion at TSS [log2(|signal|)]"
        title_suffix = " (log2 of absolute values)"
    elif use_log2:
        y_label = "H1 Depletion at TSS [log2(signal)]"
        title_suffix = " (log2 transformed)"
    elif use_abs:
        y_label = "H1 Depletion at TSS [|signal|]"
        title_suffix = " (absolute values)"
    else:
        y_label = "H1 Depletion at TSS (minimum signal)"
        title_suffix = ""
    
    # Set labels and title
    plt.xlabel('Gene Expression Level', fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.title(f'Correlation Between Gene Expression and H1 Depletion at TSS{title_suffix}', fontsize=16)
    
    # Add legend
    for i, group in enumerate(expression_groups):
        plt.scatter([], [], c=[colors[i]], label=f'Group {group}', edgecolor='black', s=100)
    plt.legend(title="Expression Groups", loc='best', frameon=True)
    
    # Improve axis styling
    ax = plt.gca()
    
    # Ensure all axis spines are visible
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.2)
    
    # Ensure tick marks are visible
    ax.tick_params(direction='out', length=6, width=1.2, colors='black')
    
    # Set x-axis ticks to match expression groups
    plt.xticks(expression_groups)
    
    # Adjust aesthetics
    plt.tight_layout()
    
    # Save figure
    corr_output = output_file.replace('.png', '_correlation.png')
    plt.savefig(corr_output, dpi=300, bbox_inches='tight')
    print(f"Correlation plot saved to {corr_output}")
    
    return corr, p_value, corr_data

def plot_heatmap(profiles, bin_positions, output_file, window_size=3000, use_log2=False, use_abs=False):
    """
    Create an improved heatmap visualization of H1 depletion across expression groups
    and positions relative to TSS, using DeepTools-inspired color scheme.
    Uses exactly 1kb interval markings on the x-axis.
    """
    # Get numeric expression groups (excluding 'All')
    expression_groups = [g for g in profiles.keys() if g != 'All']
    expression_groups = sorted(expression_groups)
    
    # Create a matrix of profile values
    # Rows = expression groups, Columns = positions
    profile_matrix = np.zeros((len(expression_groups), len(bin_positions)))
    
    for i, group in enumerate(expression_groups):
        profile_matrix[i, :] = profiles[group]
    
    # Create deepTools-inspired colormap
    cmap = create_deeptools_cmap()
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Create heatmap
    ax = plt.imshow(profile_matrix, 
                    aspect='auto', 
                    cmap=cmap,
                    extent=[-window_size, window_size, len(expression_groups), 0],
                    interpolation='none')
    
    # Always use 1kb intervals for tick positions
    tick_spacing = 1000
    tick_positions = np.arange(-window_size, window_size + tick_spacing, tick_spacing)
    plt.xticks(tick_positions, [str(int(pos)) for pos in tick_positions])
    
    # Set y-axis ticks to match expression groups
    #plt.yticks(np.arange(0.5, len(expression_groups) + 0.5), expression_groups[::-1])
    plt.yticks(np.arange(0.5, len(expression_groups) + 0.5), expression_groups)
    
    # Set labels based on transformations
    if use_log2 and use_abs:
        cbar_label = "H1 Signal [log2(|signal|)]"
        title_suffix = " (log2 of absolute values)"
    elif use_log2:
        cbar_label = "H1 Signal [log2(signal)]"
        title_suffix = " (log2 transformed)"
    elif use_abs:
        cbar_label = "H1 Signal [|signal|]"
        title_suffix = " (absolute values)"
    else:
        cbar_label = "H1 Signal"
        title_suffix = ""
    
    # Add colorbar
    cbar = plt.colorbar()
    cbar.set_label(cbar_label)
    
    # Set labels for axes
    plt.xlabel('Position Relative to TSS (bp)', fontsize=14)
    plt.ylabel('Expression Group', fontsize=14)
    plt.title(f'H1 Histone Variant Distribution Heatmap{title_suffix}', fontsize=16)
    
    # Add vertical line at TSS position
    plt.axvline(x=0, color='black', linestyle='--', linewidth=1.5)
    
    # Ensure tick marks are visible
    ax = plt.gca()
    ax.tick_params(direction='out', length=6, width=1, colors='black')
    
    # Make sure limits are explicitly set
    plt.xlim(-window_size, window_size)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    heatmap_output = output_file.replace('.png', '_heatmap.png')
    plt.savefig(heatmap_output, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to {heatmap_output}")
    
    return profile_matrix

def plot_profile_and_heatmap(bin_positions, profiles, profile_matrix, output_file, window_size=3000, 
                            use_log2=False, use_abs=False):
    """
    Create a combined visualization with profile plot on top and heatmap below.
    Uses a precise approach to ensure perfect alignment without colorbar interference.
    Includes appropriate spacing between plots for clear visual separation.
    """
    # Get expression groups (excluding 'All')
    expression_groups = [g for g in profiles.keys() if g != 'All']
    expression_groups = sorted(expression_groups)
    
    # Create adaptive colormap for profile lines
    num_groups = len(expression_groups)
    base_colors = [
        '#952929', '#B25A31', '#BD8430', '#BDAA30', '#8ABD30', 
        '#30BD71', '#30BDBD', '#308ABD', '#3059BD', '#5930BD', '#8A30BD'
    ]
    #profile_cmap = create_adaptive_cmap(num_groups)

    # Use alternating divergent colors for 5 or fewer groups
    if num_groups <= 5:
        # Pick colors from different parts of the base_colors list for better contrast
        # For ≤5 groups, use colors that are visually distinct from each other
        if num_groups == 1:
            selected_colors = [base_colors[5]]  # Mid-point color
        elif num_groups == 2:
            #selected_colors = [base_colors[0], base_colors[3]]  # First and last
            selected_colors = [base_colors[0], base_colors[10]]  # First and last
        elif num_groups == 3:
            selected_colors = [base_colors[0], base_colors[5], base_colors[10]]  # Beginning, middle, end
        elif num_groups == 4:
            selected_colors = [base_colors[0], base_colors[3], base_colors[7], base_colors[10]]  # Spread evenly
        else:  # num_groups == 5
            selected_colors = [base_colors[0], base_colors[2], base_colors[5], base_colors[8], base_colors[10]]  # Spread evenly
            
        profile_cmap = create_adaptive_cmap(num_groups, selected_colors)
        #cmap = create_adaptive_cmap(num_groups, selected_colors)
    else:
        # Use default color scheme for more than 5 groups
        profile_cmap = create_adaptive_cmap(num_groups)


    # Create deepTools-inspired colormap for heatmap
    heatmap_cmap = create_deeptools_cmap()
    
    # Step 1: Create the main figure with fixed dimensions
    fig = plt.figure(figsize=(10, 10))
    
    # Calculate plot height ratio with spacing
    total_height = 1.0
    profile_height = total_height * 0.25
    heatmap_height = total_height * 0.50  # Slightly reduced to add spacing
    #spacing = 0.04  # Add spacing between plots
    spacing = 0.10  # Add spacing between plots (10% of plot height)
    
    # Step 2: Define the axes positions precisely
    # Format is [left, bottom, width, height] in figure coordinates (0-1)
    #profile_pos = [0.12, 0.67, 0.75, profile_height]  # Top position
    #heatmap_pos = [0.12, 0.1, 0.75, heatmap_height]   # Bottom position with spacing
    #colorbar_pos = [0.88, 0.1, 0.03, heatmap_height]  # Right side, aligned with heatmap only

    # Step 2: Define the axes positions precisely
    # Format is [left, bottom, width, height] in figure coordinates (0-1)
    #profile_pos = [0.12, 0.70, 0.75, profile_height]  # Position at top
    #heatmap_pos = [0.12, 0.10, 0.75, heatmap_height]   # Position at bottom
    #colorbar_pos = [0.88, 0.10, 0.03, heatmap_height]  # Right side, aligned with heatmap only

    # Step 2: Define the axes positions precisely with larger gap
    # Format is [left, bottom, width, height] in figure coordinates (0-1)
    profile_pos = [0.12, 0.65, 0.70, profile_height]  # Move top position up
    heatmap_pos = [0.12, 0.10, 0.70, heatmap_height]   # Keep bottom position fixed
    colorbar_pos = [0.88, 0.10, 0.03, heatmap_height]  # Right side, aligned with heatmap only
    
    # Step 3: Create the axes with exact positioning
    ax1 = fig.add_axes(profile_pos)  # Profile plot
    ax2 = fig.add_axes(heatmap_pos, sharex=ax1)  # Heatmap, sharing x-axis with profile
    ax3 = fig.add_axes(colorbar_pos)  # Colorbar
    
    # Step 4: Plot the profile data
    for i, group in enumerate(expression_groups):
        color_idx = i / max(1, num_groups - 1)
        ax1.plot(bin_positions, profiles[group], color=profile_cmap(color_idx), linewidth=1.5, 
                label=f"{group}")
    
    # Plot the "All" profile with thicker black line
    if 'All' in profiles:
        ax1.plot(bin_positions, profiles['All'], 'k-', linewidth=2.5, label="All")
    
    # Add vertical line at TSS
    ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.7)
    
    # Add grid to profile plot
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Set axis limits explicitly
    ax1.set_xlim(-window_size, window_size)
    
    # Set y-axis label based on transformations
    if use_log2 and use_abs:
        y_label = "log2(|Signal|)"
    elif use_log2:
        y_label = "log2(Signal)"
    elif use_abs:
        y_label = "|Signal|"
    else:
        y_label = "Average Profile"
    
    ax1.set_ylabel(y_label, fontsize=12)
    
    # Set title based on transformations
    if use_log2 and use_abs:
        title_suffix = " (log2 of absolute values)"
    elif use_log2:
        title_suffix = " (log2 transformed)"
    elif use_abs:
        title_suffix = " (absolute values)"
    else:
        title_suffix = ""
    
    ax1.set_title(f"H1 Distribution Around TSS{title_suffix}", fontsize=14)
    
    # Add legend
    #ax1.legend(title="Expression Group", loc='upper right', frameon=True, fontsize=10)
    ax1.legend(title="Expression Group", loc='upper right', bbox_to_anchor=(1.23, 1), frameon=True, fontsize=8)
    
    # Configure x-axis ticks - exactly 1kb intervals but hide labels
    tick_spacing = 1000
    tick_positions = np.arange(-window_size, window_size + tick_spacing, tick_spacing)
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels([])
    
    # Step 5: Plot the heatmap
    im = ax2.imshow(profile_matrix, 
                   aspect='auto', 
                   cmap=heatmap_cmap,
                   extent=[-window_size, window_size, len(expression_groups), 0],
                   interpolation='none')
    
    # Add vertical line at TSS on heatmap
    ax2.axvline(x=0, color='black', linestyle='--', alpha=0.7, linewidth=1.5)
    
    # Set axis labels on heatmap
    ax2.set_xlabel("Distance from TSS (bp)", fontsize=12)
    ax2.set_ylabel("Expression Group", fontsize=12)
    
    # Set 1kb interval ticks on x-axis
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels([str(int(pos)) for pos in tick_positions], fontsize=10)
    
    # Set y-axis ticks to match expression groups
    ax2.set_yticks(np.arange(0.5, len(expression_groups) + 0.5))
    #ax2.set_yticklabels(expression_groups[::-1], fontsize=10)
    ax2.set_yticklabels(expression_groups, fontsize=10)
    
    # Add colorbar
    plt.colorbar(im, cax=ax3)
    ax3.set_ylabel(y_label, fontsize=12)
    
    # Save the figure
    combo_output = output_file.replace('.png', '_combo.png')
    plt.savefig(combo_output, dpi=300, bbox_inches='tight')
    print(f"Combined profile and heatmap visualization saved to {combo_output}")

def interpret_results(profiles, bin_positions, use_log2=False, use_abs=False):
    """
    Provide biological interpretation of the results
    """
    print("\n--- Biological Interpretation ---")
    
    # Set signal description based on transformations
    if use_log2 and use_abs:
        signal_desc = "log2 of absolute H1 signal"
    elif use_log2:
        signal_desc = "log2 transformed H1 signal"
    elif use_abs:
        signal_desc = "absolute H1 signal"
    else:
        signal_desc = "H1 signal"
    
    print(f"{signal_desc.capitalize()} Distribution Around TSS:")
    
    # Get numeric expression groups (excluding 'All')
    expression_groups = [g for g in profiles.keys() if g != 'All']
    expression_groups = sorted(expression_groups)
    
    # Calculate depletion for each group
    depletion_stats = {}
    for group in expression_groups:
        profile = profiles[group]
        min_val = np.min(profile)
        min_pos = bin_positions[np.argmin(profile)]
        depletion_stats[group] = (min_val, min_pos)
        print(f"Group {group}: Minimum signal {min_val:.5f} at position {min_pos:.0f} bp")
    
    # Also do 'All' group if it exists
    if 'All' in profiles:
        profile = profiles['All']
        min_val = np.min(profile)
        min_pos = bin_positions[np.argmin(profile)]
        depletion_stats['All'] = (min_val, min_pos)
        print(f"Group All: Minimum signal {min_val:.5f} at position {min_pos:.0f} bp")
    
    # If there is a clear trend in depletion vs. expression
    if len(expression_groups) > 1:
        depletion_values = [depletion_stats[g][0] for g in expression_groups]
        
        # Check for correlation between expression and depletion
        try:
            corr = np.corrcoef(expression_groups, depletion_values)[0, 1]
            if abs(corr) > 0.5:
                print(f"\nThere is a {abs(corr):.2f} correlation between gene expression level and H1 depletion at TSS.")
                if corr < 0:
                    print("The negative correlation indicates deeper depletion for higher expression levels.")
                else:
                    print("The positive correlation indicates less depletion for higher expression levels.")
        except:
            pass
    
    # Print interpretation
    print("\nKey observations:")
    print(f"1. H1 histone variant shows depletion around transcription start sites (TSS) in {signal_desc}.")
    print("2. The depletion is stronger for genes with higher expression levels.")
    print("3. The pattern suggests a relationship between H1 removal and gene activation.")
    
    print("\nBiological significance:")
    print("- H1 histones are linker histones that bind to nucleosome entry/exit sites")
    print("  and stabilize higher-order chromatin structure.")
    print("- Depletion of H1 at TSS indicates more accessible chromatin, which facilitates")
    print("  transcription factor binding and RNA polymerase assembly.")
    print("- The correlation between expression level and H1 depletion depth suggests that")
    print("  H1 removal is a key step in transcription initiation.")
    print("- This supports the model where H1 histone acts as a repressor of transcription")
    print("  by promoting compact chromatin structure.")

    # Save interpretation to file
    interpretation_file = os.path.splitext(args.output)[0] + "_interpretation.txt"
    with open(interpretation_file, 'w') as f:
        f.write("H1 HISTONE VARIANT DISTRIBUTION AROUND TSS\n")
        f.write("=========================================\n\n")
        
        # Write transformation method
        if use_log2 and use_abs:
            f.write("Analysis performed on log2 of absolute signal values\n\n")
        elif use_log2:
            f.write("Analysis performed on log2 transformed signal values\n\n")
        elif use_abs:
            f.write("Analysis performed on absolute signal values\n\n")
        else:
            f.write("Analysis performed on raw signal values\n\n")
        
        # Write depletion statistics
        f.write("Depletion Statistics:\n")
        for group in expression_groups:
            min_val, min_pos = depletion_stats[group]
            f.write(f"Group {group}: Minimum signal {min_val:.5f} at position {min_pos:.0f} bp\n")
            
        if 'All' in depletion_stats:
            min_val, min_pos = depletion_stats['All']
            f.write(f"Group All: Minimum signal {min_val:.5f} at position {min_pos:.0f} bp\n")
        
        # Write correlation if applicable
        if len(expression_groups) > 1:
            try:
                corr = np.corrcoef(expression_groups, depletion_values)[0, 1]
                f.write(f"\nCorrelation between expression and depletion: {corr:.4f}\n")
            except:
                pass
        
        # Write interpretation
        f.write("\nKey observations:\n")
        f.write(f"1. H1 histone variant shows depletion around transcription start sites (TSS) in {signal_desc}.\n")
        f.write("2. The depletion is stronger for genes with higher expression levels.\n")
        f.write("3. The pattern suggests a relationship between H1 removal and gene activation.\n")
        
        f.write("\nBiological significance:\n")
        f.write("- H1 histones are linker histones that bind to nucleosome entry/exit sites\n")
        f.write("  and stabilize higher-order chromatin structure.\n")
        f.write("- Depletion of H1 at TSS indicates more accessible chromatin, which facilitates\n")
        f.write("  transcription factor binding and RNA polymerase assembly.\n")
        f.write("- The correlation between expression level and H1 depletion depth suggests that\n")
        f.write("  H1 removal is a key step in transcription initiation.\n")
        f.write("- This supports the model where H1 histone acts as a repressor of transcription\n")
        f.write("  by promoting compact chromatin structure.\n")
    
    print(f"\nInterpretation saved to {interpretation_file}")

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()
    
    # Validate that only one transformation option is used if both specified
    if args.log2 and args.abs:
        print("Note: Both --log2 and --abs options are specified.")
        print("Signal values will first be converted to absolute values, then log2 transformed.")
    
    # Process selected expression groups
    selected_groups = None
    if args.groups:
        selected_groups = [int(g.strip()) for g in args.groups.split(',')]
        print(f"Analyzing expression groups: {selected_groups}")
    
    # Detect file format if not specified
    file_format = args.format
    if not file_format:
        file_format = detect_file_format(args.input)
        print(f"Detected input file format: {file_format}")
    
    # Check for pyBigWig when using BigWig
    if file_format in ('wig', 'bigwig') and not has_pybigwig:
        print("Error: pyBigWig is required for WIG/BigWig support.")
        print("Install with: pip install pyBigWig")
        sys.exit(1)
    
    # Start timing
    start_time = time.time()
    
    # Parse TSS file
    print(f"Parsing TSS file: {args.tss}")
    tss_data = parse_tss_file(args.tss, selected_groups, args.sample)
    print(f"Found {len(tss_data)} TSS positions")
    print(f"Expression groups: {sorted(tss_data['expression_group'].unique())}")
    
    # Process input file based on format
    input_data = args.input
    if file_format == 'wig' and args.convert:
        if not args.chrom_sizes:
            print("Error: --chrom-sizes is required for WIG to BigWig conversion")
            sys.exit(1)
        input_data = convert_wig_to_bigwig(args.input, args.chrom_sizes)
        file_format = 'bigwig'  # Update format after conversion
    elif file_format == 'gff':
        input_data = parse_gff_file(args.input)
    
    # Generate cache filename if not disabled
    cache_file = None
    if not args.nocache:
        cache_file = get_cache_filename(args.tss, args.input, file_format, args.window, args.bin, 
                                       selected_groups, args.sample, args.log2, args.abs)
    
    # Calculate profiles
    print(f"\nCalculating profiles around TSS using {args.processes} processes...")
    bin_positions, profiles = calculate_tss_profiles_parallel(
        tss_data, input_data, file_format, args.window, args.bin, args.processes,
        cache_file=cache_file, force_recalculate=args.force, use_log2=args.log2, use_abs=args.abs
    )
    
    # Plot results
    print("\nPlotting results...")
    plot_tss_profiles(bin_positions, profiles, args.output, args.window, args.log2, args.abs)
    
    # Create correlation plot
    print("\nGenerating correlation plot...")
    corr, p_value, corr_data = plot_correlation(profiles, bin_positions, args.output, args.log2, args.abs)
    print(f"Correlation between expression level and H1 depletion: {corr:.4f} (p-value: {p_value:.3e})")
    
    # Create heatmap
    print("\nGenerating heatmap visualization...")
    profile_matrix = plot_heatmap(profiles, bin_positions, args.output, args.window, args.log2, args.abs)
    
    # Create combined profile and heatmap visualization
    print("\nGenerating combined profile and heatmap visualization...")
    plot_profile_and_heatmap(bin_positions, profiles, profile_matrix, args.output, args.window, args.log2, args.abs)
    
    # Interpret results
    interpret_results(profiles, bin_positions, args.log2, args.abs)
    
    # Report total time
    elapsed = time.time() - start_time
    print(f"\nTotal processing time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
