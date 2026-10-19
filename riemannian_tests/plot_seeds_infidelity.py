#!/usr/bin/env python3
"""
Plot infidelity convergence from all seed_* directories.

This script automatically finds all seed_* directories, reads optim_history.dat
from each, and plots the infidelity convergence on a log scale.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
import argparse

def find_seed_directories(parent_dir='.', basename='seed'):
    """
    Find all directories matching the basename pattern in the parent directory.

    Parameters
    ----------
    parent_dir : str
        Parent directory to search (default: current directory)
    basename : str
        Basename pattern for seed directories (default: 'seed')
        Looks for directories matching: {basename}_*

    Returns
    -------
    list of str
        Sorted list of seed directory paths
    """
    pattern = os.path.join(parent_dir, f'{basename}_*')
    seed_dirs = sorted(glob.glob(pattern))

    # Filter to only actual directories
    seed_dirs = [d for d in seed_dirs if os.path.isdir(d)]

    return seed_dirs

def read_infidelity(seed_dir):
    """
    Read infidelity data from optim_history.dat in seed directory.

    Parameters
    ----------
    seed_dir : str
        Path to seed directory

    Returns
    -------
    tuple (iteration, infidelity, seed_num) or None if file not found/readable
    """
    history_file = os.path.join(seed_dir, 'optim_history.dat')

    if not os.path.exists(history_file):
        print(f"Warning: {history_file} not found, skipping")
        return None

    try:
        data = np.loadtxt(history_file, skiprows=1)

        if data.ndim == 1:
            # Single iteration
            data = data.reshape(1, -1)

        iteration = data[:, 0].astype(int)
        infidelity = data[:, 5]  # Column 6 (1-indexed) = index 5 (0-indexed)

        # Extract seed number from directory name
        # Pattern: {basename}_XXXX -> extract XXXX as integer
        seed_name = os.path.basename(seed_dir)
        # Get the part after the last underscore
        seed_num = int(seed_name.split('_')[-1])

        return (iteration, infidelity, seed_num)

    except Exception as e:
        print(f"Warning: Error reading {history_file}: {e}")
        return None

def plot_seeds_infidelity(parent_dir='.', basename='seed', title_suffix='', show_individual=True, show_statistics=True):
    """
    Plot infidelity convergence from all seed directories.

    Statistics are computed by interpolating on an infidelity grid:
    - For each infidelity level, find the iteration number where it's reached
    - Compute mean and std dev of iteration numbers across seeds

    Parameters
    ----------
    parent_dir : str
        Parent directory containing seed directories
    basename : str
        Basename pattern for seed directories (default: 'seed')
    title_suffix : str
        Optional suffix to add to title after directory name (default: '')
    show_individual : bool
        If True, plot individual seed traces (default: True)
    show_statistics : bool
        If True, plot mean and std bands (default: True)
    """
    # Find all seed directories
    seed_dirs = find_seed_directories(parent_dir, basename)

    if len(seed_dirs) == 0:
        print(f"Error: No {basename}_* directories found in {parent_dir}")
        return

    print(f"Found {len(seed_dirs)} seed directories")

    # Read data from all seeds
    all_data = []
    for seed_dir in seed_dirs:
        result = read_infidelity(seed_dir)
        if result is not None:
            all_data.append(result)

    if len(all_data) == 0:
        print("Error: No valid data found in any seed directory")
        return

    print(f"Successfully read data from {len(all_data)} seeds")

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Plot individual seeds
    if show_individual:
        for iteration, infidelity, seed_num in all_data:
            ax.plot(iteration, infidelity,
                   linestyle='-', linewidth=0.8, alpha=0.3,
                   color='C0', label='_nolegend_')

        # Add one labeled line for the legend
        ax.plot([], [], linestyle='-', linewidth=0.8, alpha=0.3,
               color='C0', label=f'Individual seeds (n={len(all_data)})')

    # Compute and plot statistics on interpolated grid
    if show_statistics and len(all_data) > 1:
        # Find global infidelity range across all seeds
        all_infidelities = np.concatenate([data[1] for data in all_data])
        min_infidelity_global = np.min(all_infidelities)
        max_infidelity_global = np.max(all_infidelities)

        # Create infidelity grid (logarithmically spaced)
        n_grid = 100
        infidelity_grid = np.logspace(np.log10(min_infidelity_global),
                                      np.log10(max_infidelity_global),
                                      n_grid)

        # For each seed, interpolate to find iteration at each infidelity level
        iteration_at_infidelity = []

        for iteration, infidelity, seed_num in all_data:
            # Infidelity should be monotonically decreasing, but might have small increases
            # Use only the decreasing envelope to avoid interpolation issues

            # Create monotonic decreasing infidelity by taking cumulative minimum
            infidelity_mono = np.minimum.accumulate(infidelity)

            # Interpolate in log space (more stable for log-scale data)
            # Need to reverse arrays since np.interp requires x to be increasing
            iter_interp = np.interp(
                np.log10(infidelity_grid),           # x values to interpolate at
                np.log10(infidelity_mono)[::-1],     # x data (reversed - increasing)
                iteration[::-1],                      # y data (reversed)
                left=np.nan, right=np.nan             # NaN for out-of-range values
            )

            iteration_at_infidelity.append(iter_interp)

        # Convert to array: shape (n_seeds, n_grid)
        iteration_at_infidelity = np.array(iteration_at_infidelity)

        # Compute statistics (ignoring NaN values)
        with np.errstate(invalid='ignore'):  # Suppress warnings for all-NaN slices
            mean_iterations = np.nanmean(iteration_at_infidelity, axis=0)
            std_iterations = np.nanstd(iteration_at_infidelity, axis=0)
            min_iterations = np.nanmin(iteration_at_infidelity, axis=0)
            max_iterations = np.nanmax(iteration_at_infidelity, axis=0)

        # Only plot where we have valid data from ALL seeds
        n_valid = np.sum(~np.isnan(iteration_at_infidelity), axis=0)
        valid_mask = (n_valid == len(all_data))

        if np.any(valid_mask):
            infidelity_valid = infidelity_grid[valid_mask]
            mean_iter_valid = mean_iterations[valid_mask]
            std_iter_valid = std_iterations[valid_mask]
            min_iter_valid = min_iterations[valid_mask]
            max_iter_valid = max_iterations[valid_mask]

            # Plot mean
            ax.plot(mean_iter_valid, infidelity_valid,
                   linestyle='-', linewidth=2.5, color='red',
                   marker='o', markersize=5, markevery=max(1, len(mean_iter_valid)//20),
                   label='Mean', zorder=10)

            # Plot standard deviation band
            ax.fill_betweenx(infidelity_valid,
                           mean_iter_valid - std_iter_valid,
                           mean_iter_valid + std_iter_valid,
                           alpha=0.2, color='red', label='±1 std dev')

            # Also show min and max
            ax.plot(min_iter_valid, infidelity_valid,
                   linestyle='--', linewidth=1.5, color='green',
                   alpha=0.6, label='Min/Max', zorder=5)
            ax.plot(max_iter_valid, infidelity_valid,
                   linestyle='--', linewidth=1.5, color='green',
                   alpha=0.6, label='_nolegend_', zorder=5)

            # Set y-axis limits: start where all seeds have data, end at 1.0
            ax.set_ylim([np.min(infidelity_valid), 1.0])

    # Format plot
    # Get directory name for title
    dir_name = os.path.basename(os.path.abspath(parent_dir))

    # Build title with optional suffix
    if title_suffix:
        title = f'Infidelity Convergence: {dir_name} {title_suffix} (n={len(all_data)} seeds)'
    else:
        title = f'Infidelity Convergence: {dir_name} (n={len(all_data)} seeds)'

    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Infidelity', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')
    ax.set_yscale('log')

    # Set y-axis limits and ticks to ensure 10^-6 and 1.0 are visible
    from matplotlib.ticker import LogLocator, LogFormatter
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
    ax.yaxis.set_major_formatter(LogFormatter(base=10.0, labelOnlyBase=False))

    # Ensure y-limits include 1e-6 and 1.0
    current_ylim = ax.get_ylim()
    new_ymin = min(current_ylim[0], 1e-6)
    new_ymax = max(current_ylim[1], 1.0)
    ax.set_ylim([new_ymin, new_ymax])

    plt.tight_layout()

    # Save the figure
    n_seeds = len(all_data)
    if title_suffix:
        # Create filename-safe version of suffix (replace spaces and special chars)
        suffix_safe = title_suffix.replace(' ', '_').replace('/', '_').replace('\\', '_')
        output_filename = os.path.join(parent_dir, f'seeds_infidelity_convergence_{suffix_safe}_n{n_seeds}.png')
    else:
        output_filename = os.path.join(parent_dir, f'seeds_infidelity_convergence_n{n_seeds}.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_filename}")

    # Print summary statistics
    if show_statistics and len(all_data) > 1:
        print("\nSummary Statistics (Final Infidelity):")
        final_infidelities = [data[1][-1] for data in all_data]
        print(f"  Mean:   {np.mean(final_infidelities):.6e}")
        print(f"  Std:    {np.std(final_infidelities):.6e}")
        print(f"  Min:    {np.min(final_infidelities):.6e}")
        print(f"  Max:    {np.max(final_infidelities):.6e}")
        print(f"  Median: {np.median(final_infidelities):.6e}")

    # Show the plot
    print("\nDisplaying plot... Close the window to exit.")
    plt.show(block=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plot infidelity convergence from multiple seed directories.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                                # Use current directory, default basename
  %(prog)s qft2                                           # Specify parent directory
  %(prog)s qft2 --suffix "MINRES+LR"                      # Add title suffix
  %(prog)s . --basename "config_QFT2_dual_seed"           # Custom seed directory basename
  %(prog)s . --basename "config_QFT2_lbfgs_seed" --suffix "L-BFGS"

Statistics:
  Mean and standard deviation are computed by interpolating on an infidelity grid.
  For each infidelity level, the iteration number where it's reached is found for
  each seed, then statistics are computed across seeds. Only infidelity levels where
  ALL seeds have data are shown.

Output:
  - Plot displayed interactively
  - PNG file saved in parent directory with naming:
      seeds_infidelity_convergence_[suffix_]n<N>.png
  - Summary statistics printed to console
        """
    )

    parser.add_argument(
        'directory',
        nargs='?',
        default='.',
        help='Parent directory containing seed subdirectories (default: current directory)'
    )

    parser.add_argument(
        '--basename',
        type=str,
        default='seed',
        help='Basename pattern for seed directories (default: "seed", looks for seed_*)'
    )

    parser.add_argument(
        '--suffix',
        type=str,
        default='',
        help='Optional suffix to add to title and filename (e.g., "MINRES+Low-Rank")'
    )

    parser.add_argument(
        '--no-individual',
        action='store_true',
        help='Hide individual seed traces (only show statistics)'
    )

    parser.add_argument(
        '--no-statistics',
        action='store_true',
        help='Hide statistics (only show individual traces)'
    )

    args = parser.parse_args()

    # Check if parent directory exists
    if not os.path.isdir(args.directory):
        print(f"Error: Directory '{args.directory}' does not exist")
        sys.exit(1)

    plot_seeds_infidelity(
        args.directory,
        args.basename,
        args.suffix,
        show_individual=not args.no_individual,
        show_statistics=not args.no_statistics
    )
