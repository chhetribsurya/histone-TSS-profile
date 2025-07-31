#!/usr/bin/env python3
"""
Smoothing filter for H1 histone variant distribution around TSS.
This script adds smoothing to existing TSS profile data and creates
visualizations with original and smoothed profiles side by side.

Usage:
    python smooth_tss_profiles.py --input profile_data.csv --output smoothed_profiles.png

Input CSV should have columns: 'Position', 'Group_1', 'Group_2', ..., 'Group_All'
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
import seaborn as sns

# Set matplotlib parameters
# plt.rcParams['font.size'] = 10
# plt.rcParams['figure.figsize'] = (5.6, 3.8)

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Apply smoothing to TSS profiles')
    parser.add_argument('--input', required=True, help='Input CSV file with profile data')
    parser.add_argument('--output-dir', required=True, help='Directory to save all output files')
    parser.add_argument('--prefix', default='smoothed', help='Prefix for output files')
    parser.add_argument('--window', type=int, default=3000, help='Window size around TSS (bp)')
    
    # Smoothing options
    parser.add_argument('--method', choices=['savgol', 'gaussian', 'moving_avg'], 
                     default='savgol', help='Smoothing method to use')
    parser.add_argument('--smooth-window', type=int, default=21, 
                     help='Window size for smoothing (Savitzky-Golay or moving average)')
    parser.add_argument('--smooth-poly', type=int, default=3, 
                     help='Polynomial order for Savitzky-Golay filter')
    parser.add_argument('--smooth-sigma', type=float, default=2, 
                     help='Sigma (standard deviation) for Gaussian filter')
    
    # Visualization options
    parser.add_argument('--log2', action='store_true', help='Data is log2 transformed')
    parser.add_argument('--abs', action='store_true', help='Data is absolute values')
    parser.add_argument('--no-density', action='store_true', help='Skip density plot generation')
    
    return parser.parse_args()

def load_profile_data(csv_file):
    """Load profile data from CSV file"""
    try:
        df = pd.read_csv(csv_file)
        
        # Extract bin positions
        if 'Position' in df.columns:
            bin_positions = df['Position'].values
        else:
            # If no Position column, assume first column is positions
            bin_positions = df.iloc[:, 0].values
            
        # Extract profiles for each group
        profiles = {}
        for col in df.columns:
            if col.startswith('Group_') or col == 'All':
                group = col.replace('Group_', '')
                try:
                    # Try to convert group to int if it's numeric
                    group = int(group)
                except ValueError:
                    pass
                
                profiles[group] = df[col].values
        
        return bin_positions, profiles
    
    except Exception as e:
        print(f"Error loading profile data: {e}")
        sys.exit(1)

def apply_smoothing(profiles, bin_positions, method='savgol', window=21, polyorder=3, sigma=2):
    """Apply smoothing to profile data using various methods"""
    smoothed_profiles = {}
    
    # For each expression group
    for group in profiles:
        profile = profiles[group]
        
        # Apply specified smoothing method
        if method == 'savgol':
            # Savitzky-Golay filter (requires odd window size)
            if window % 2 == 0:
                window += 1
            
            # Ensure polyorder < window
            if polyorder >= window:
                polyorder = window - 1
                
            smoothed_profile = savgol_filter(profile, window, polyorder)
        
        elif method == 'gaussian':
            # Gaussian filter
            smoothed_profile = gaussian_filter1d(profile, sigma)
        
        elif method == 'moving_avg':
            # Simple moving average
            kernel = np.ones(window) / window
            smoothed_profile = np.convolve(profile, kernel, mode='same')
            
            # Handle edges
            half_window = window // 2
            for i in range(half_window):
                # Left edge
                edge_kernel = np.ones(i + half_window + 1) / (i + half_window + 1)
                smoothed_profile[i] = np.sum(profile[:i + half_window + 1] * edge_kernel)
                
                # Right edge
                j = len(profile) - i - 1
                edge_kernel = np.ones(i + half_window + 1) / (i + half_window + 1)
                smoothed_profile[j] = np.sum(profile[j - half_window:] * edge_kernel)
        
        else:
            raise ValueError(f"Unknown smoothing method: {method}")
        
        smoothed_profiles[group] = smoothed_profile
    
    return smoothed_profiles

def plot_smoothed_tss_profiles(bin_positions, profiles, smoothed_profiles, output_file, 
                              window_size=3000, method='savgol', params=None, use_log2=False, use_abs=False):
    """
    Create a visualization with both original and smoothed profiles side by side.
    Saves the output file to the specified path.
    """
    plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 8,
    'figure.titlesize': 12,
    'figure.figsize': (5.6, 3.8)
    })
    
    # Get expression groups (excluding 'All')
    expression_groups = [g for g in profiles.keys() if g != 'All']
    expression_groups = sorted(expression_groups)
    
    # Define colors for highly and lowly expressed enhancers
    base_colors = ['orangered', 'gray']  # Highly expressed: orangered, Lowly expressed: gray
    
    # Get number of groups
    num_groups = len(expression_groups)
    
    # Create color mapping
    group_colors = {}
    for i, group in enumerate(expression_groups):
        # Use orangered for highly expressed and gray for lowly expressed
        group_colors[group] = base_colors[0] if i >= num_groups/2 else base_colors[1]
    
    # Create figure with side-by-side panels for original and smoothed
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1], width_ratios=[1, 1])
    
    # Left panel: Original profiles
    ax1 = plt.subplot(gs[0, 0])
    
    # Plot each expression group using the assigned colors
    for group in expression_groups:
        ax1.plot(bin_positions, profiles[group], color=group_colors[group], linewidth=1, 
                 label=f"{group}")
    
    # Plot the "All" profile with thicker black line
    if 'All' in profiles:
        ax1.plot(bin_positions, profiles['All'], 'k-', linewidth=2, label="All")
    
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
    ax1.set_title("Original Profile")
    
    # Right panel: Smoothed profiles
    ax2 = plt.subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    
    # Plot each expression group using the assigned colors
    for group in expression_groups:
        ax2.plot(bin_positions, smoothed_profiles[group], color=group_colors[group], linewidth=1.5, 
                 label=f"{group}")
    
    # Plot the "All" profile with thicker black line
    if 'All' in smoothed_profiles:
        ax2.plot(bin_positions, smoothed_profiles['All'], 'k-', linewidth=2.5, label="All")
    
    # Add vertical line at TSS
    ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.7)
    
    # Add grid
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    # Set smoothing method information
    if method == 'savgol':
        method_str = f"Savitzky-Golay (window={params['window']}, order={params['polyorder']})"
    elif method == 'gaussian':
        method_str = f"Gaussian (sigma={params['sigma']})"
    elif method == 'moving_avg':
        method_str = f"Moving Average (window={params['window']})"
    else:
        method_str = method
    
    ax2.set_title("H1.0 Distribution")
    
    # Set overall title based on transformations
    if use_log2 and use_abs:
        title_suffix = " (log2 of absolute values)"
    elif use_log2:
        title_suffix = " (log2 transformed)"
    elif use_abs:
        title_suffix = " (absolute values)"
    else:
        title_suffix = ""
    
    fig.suptitle(f"Distribution Around TSS{title_suffix}", fontsize=16)
    
    # Bottom panels: Color legend
    ax3 = plt.subplot(gs[1, 0:2])
    
    # Create a horizontal colorbar for expression groups
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    
    # Plot colorbar with the appropriate colormap
    ax3.imshow(gradient, aspect='auto', cmap=LinearSegmentedColormap.from_list('expression_cmap', base_colors, N=2))
    
    # Configure legend axis
    ax3.set_yticks([])
    
    # Set x-axis ticks to match expression groups
    ax3.set_xticks(np.linspace(0, 255, num_groups))
    
    ax3.set_xticklabels(expression_groups)
    
    # Add "All" marker
    ax3.text(265, 0.5, "All", ha='left', va='center', fontweight='bold')
    ax3.plot([260], [0.5], 'ko', markersize=8)
    
    # Add labels for x-axis and colorbar
    ax1.set_xlabel("Relative Distance (bp)")
    ax2.set_xlabel("Relative Distance (bp)")
    ax3.set_xlabel("Gene expression")
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    # Also save as SVG
    svg_file = output_file.replace('.png', '.svg')
    plt.savefig(svg_file, bbox_inches='tight')
    print(f"Combined figure saved to {output_file} and {svg_file}")
    
    # CSV file path will be managed by the main function
    # Get the directory and use it for the CSV path
    output_dir = os.path.dirname(output_file)
    csv_file = os.path.join(output_dir, os.path.basename(output_file).replace('.png', '_data.csv').replace('.svg', '_data.csv'))
    data_for_csv = {'Position': bin_positions}
    
    # Add profiles to CSV data
    for group in expression_groups:
        data_for_csv[f'Original_Group_{group}'] = profiles[group]
        data_for_csv[f'Smoothed_Group_{group}'] = smoothed_profiles[group]
    
    if 'All' in profiles:
        data_for_csv['Original_Group_All'] = profiles['All']
        data_for_csv['Smoothed_Group_All'] = smoothed_profiles['All']
    
    pd.DataFrame(data_for_csv).to_csv(csv_file, index=False)
    print(f"Profile data saved to {csv_file}")

def plot_smoothed_kde(bin_positions, profiles, smoothed_profiles, output_file, 
                     window_size=3000, method='gaussian', use_log2=False, use_abs=False):
    """
    Create KDE-like visualization of H1 distribution around TSS.
    This provides a density plot representation of the profile data.
    Saves the output file to the specified path.
    """
    plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 8,
    'figure.titlesize': 12,
    'figure.figsize': (5.6, 3.8)
    })

    # Get expression groups (excluding 'All')
    expression_groups = [g for g in profiles.keys() if g != 'All']
    expression_groups = sorted(expression_groups)
    
    # Get number of groups
    num_groups = len(expression_groups)
    
    # Define colors for highly and lowly expressed enhancers
    base_colors = ['orangered', 'gray']  # Highly expressed: orangered, Lowly expressed: gray
    
    # Create figure
    plt.figure()
    
    # Calculate appropriate alpha for overlapping contours
    alpha = max(0.1, min(0.7, 1.0 / len(expression_groups)))
    
    # Create color mapping
    group_colors = {}
    for i, group in enumerate(expression_groups):
        # Use orangered for highly expressed and gray for lowly expressed
        group_colors[group] = base_colors[0] if i >= num_groups/2 else base_colors[1]
    
    # Plot smoothed densities using the assigned colors
    for group in expression_groups:
        # Get profile and color for this group
        profile = smoothed_profiles[group]
        color = group_colors[group]
        
        # Create density-like representation
        plt.fill_between(bin_positions, 0, profile, 
                        color=color, alpha=alpha, label=f"Group {group}")
        
        # Add contour line
        plt.plot(bin_positions, profile, color=color, linewidth=1.5)
    
    # Plot the "All" profile with thicker black line
    if 'All' in smoothed_profiles:
        plt.plot(bin_positions, smoothed_profiles['All'], 'k-', linewidth=2.5, label="All")
    
    # Add vertical line at TSS
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.7)
    
    # Add grid
    plt.grid(True, alpha=0.3, linestyle='--')
    
    # Set axis limits
    plt.xlim(-window_size, window_size)
    
    # Set y-axis label based on transformations
    if use_log2 and use_abs:
        y_label = "Density of log2(|Signal|)"
    elif use_log2:
        y_label = "Density of log2(Signal)"
    elif use_abs:
        y_label = "Density of |Signal|"
    else:
        y_label = "Signal Density"
    
    plt.ylabel(y_label)
    plt.xlabel("Relative Distance (bp)")
    
    # Set title based on transformations
    if use_log2 and use_abs:
        title_suffix = " (log2 of absolute values)"
    elif use_log2:
        title_suffix = " (log2 transformed)"
    elif use_abs:
        title_suffix = " (absolute values)"
    else:
        title_suffix = ""
    
    plt.title(f"H1 Distribution Density Around TSS{title_suffix}")
    
    # Add legend
    plt.legend(title="Expression Groups", loc='upper right')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    # Also save as SVG
    svg_file = output_file.replace('.png', '.svg')
    plt.savefig(svg_file, bbox_inches='tight')
    print(f"Density plot saved to {output_file} and {svg_file}")

def plot_standalone_smoothed_profiles(bin_positions, smoothed_profiles, output_file, 
                                    window_size=3000, method='savgol', params=None, use_log2=False, use_abs=False):
    """
    Create a standalone visualization of smoothed profiles without legends.
    Uses plt.rcParams.update() for clean formatting.
    """
    # Set font sizes and figure size for standalone plot
    plt.rcParams.update({
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.figsize': (5.6, 3.8)
    })

    # Get expression groups (excluding 'All')
    expression_groups = [g for g in smoothed_profiles.keys() if g != 'All']
    expression_groups = sorted(expression_groups)
    
    # Define colors for highly and lowly expressed enhancers
    base_colors = ['orangered', 'gray']  # Highly expressed: orangered, Lowly expressed: gray
    
    # Get number of groups
    num_groups = len(expression_groups)
    
    # Create color mapping
    group_colors = {}
    for i, group in enumerate(expression_groups):
        # Use orangered for highly expressed and gray for lowly expressed
        group_colors[group] = base_colors[0] if i >= num_groups/2 else base_colors[1]
    
    # Create figure
    plt.figure()
    
    # Plot each expression group using the assigned colors
    for group in expression_groups:
        plt.plot(bin_positions, smoothed_profiles[group], color=group_colors[group], linewidth=1.5)
    
    # Plot the "All" profile with thicker black line
    if 'All' in smoothed_profiles:
        plt.plot(bin_positions, smoothed_profiles['All'], 'k-', linewidth=2.5)
    
    # Add vertical line at TSS
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.7)
    
    # Add grid
    plt.grid(True, alpha=0.3, linestyle='--')
    
    # Set axis limits and labels
    plt.xlim(-window_size, window_size)
    
    # Set y-axis label based on transformations
    if use_log2 and use_abs:
        y_label = "log2(|Signal|)"
    elif use_log2:
        y_label = "log2(Signal)"
    elif use_abs:
        y_label = "|Signal|"
    else:
        y_label = "Average Profile"
    
    plt.ylabel(y_label)
    plt.xlabel("Relative Distance (bp)")
    
    # Set title based on transformations
    if use_log2 and use_abs:
        title_suffix = " (log2 of absolute values)"
    elif use_log2:
        title_suffix = " (log2 transformed)"
    elif use_abs:
        title_suffix = " (absolute values)"
    else:
        title_suffix = ""
    
    plt.title(f"H1 Distribution Around TSS{title_suffix}")
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    # Also save as SVG
    svg_file = output_file.replace('.png', '.svg')
    plt.savefig(svg_file, bbox_inches='tight')
    print(f"Standalone smoothed plot saved to {output_file} and {svg_file}")

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"Created output directory: {args.output_dir}")
    
    # Load profile data
    print(f"Loading profile data from {args.input}...")
    bin_positions, profiles = load_profile_data(args.input)
    print(f"Loaded profiles for {len(profiles)} expression groups")
    
    # Save a copy of the input file to the output directory
    input_filename = os.path.basename(args.input)
    input_copy_path = os.path.join(args.output_dir, input_filename)
    import shutil
    shutil.copy2(args.input, input_copy_path)
    print(f"Saved copy of input file to: {input_copy_path}")
    
    # Prepare smoothing parameters
    smoothing_params = {
        'window': args.smooth_window,
        'polyorder': args.smooth_poly, 
        'sigma': args.smooth_sigma
    }
    
    # Apply smoothing
    print(f"Applying {args.method} smoothing to profiles...")
    smoothed_profiles = apply_smoothing(
        profiles, bin_positions, method=args.method,
        window=args.smooth_window, polyorder=args.smooth_poly, sigma=args.smooth_sigma
    )
    
    # Generate output file paths
    comparison_output = os.path.join(args.output_dir, f"{args.prefix}_comparison.png")
    density_output = os.path.join(args.output_dir, f"{args.prefix}_density.png")
    standalone_output = os.path.join(args.output_dir, f"{args.prefix}_standalone.svg")
    csv_output = os.path.join(args.output_dir, f"{args.prefix}_profiles.csv")
    
    # Plot original vs smoothed profiles
    print("Generating comparison plot...")
    plot_smoothed_tss_profiles(
        bin_positions, profiles, smoothed_profiles, comparison_output,
        args.window, args.method, smoothing_params, args.log2, args.abs
    )
    
    # Plot KDE-style density visualization
    if not args.no_density:
        print("Generating density visualization...")
        plot_smoothed_kde(
            bin_positions, profiles, smoothed_profiles, density_output,
            args.window, args.method, args.log2, args.abs
        )
    
    # Plot standalone smoothed profiles (no legend)
    print("Generating standalone smoothed plot...")
    plot_standalone_smoothed_profiles(
        bin_positions, smoothed_profiles, standalone_output,
        args.window, args.method, smoothing_params, args.log2, args.abs
    )
    
    print(f"Outputs saved to directory: {args.output_dir}")
    print("Smoothing and visualization complete.")

