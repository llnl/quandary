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
    tuple (iteration, infidelity, seed_num, x0_norm) or None if file not found/readable
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

        # Extract initial control vector norm (column 13 = index 12)
        x0_norm = data[0, 12] if data.shape[1] >= 13 else None

        # Extract seed number from directory name
        # Pattern: {basename}_XXXX -> extract XXXX as integer
        seed_name = os.path.basename(seed_dir)
        # Get the part after the last underscore
        seed_num = int(seed_name.split('_')[-1])

        return (iteration, infidelity, seed_num, x0_norm)

    except Exception as e:
        print(f"Warning: Error reading {history_file}: {e}")
        return None

def plot_seeds_infidelity(parent_dir='.', basenames=None, title_suffix='', show_individual=True, show_statistics=True):
    """
    Plot infidelity convergence from all seed directories.

    Statistics are computed by interpolating on an infidelity grid:
    - For each infidelity level, find the iteration number where it's reached
    - Compute mean and std dev of iteration numbers across seeds

    Parameters
    ----------
    parent_dir : str
        Parent directory containing seed directories
    basenames : list of str
        List of basename patterns for seed directories (default: ['seed'])
    title_suffix : str
        Optional suffix to add to title after directory name (default: '')
    show_individual : bool
        If True, plot individual seed traces (default: True)
    show_statistics : bool
        If True, plot mean and std bands (default: True)
    """
    if basenames is None:
        basenames = ['seed']

    # Ensure basenames is a list
    if isinstance(basenames, str):
        basenames = [basenames]

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Color cycle for different basenames
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

    # Process each basename separately
    for basename_idx, basename in enumerate(basenames):
        color = colors[basename_idx % len(colors)]

        # Find all seed directories for this basename
        seed_dirs = find_seed_directories(parent_dir, basename)

        if len(seed_dirs) == 0:
            print(f"Warning: No {basename}_* directories found in {parent_dir}")
            continue

        print(f"Found {len(seed_dirs)} directories for basename '{basename}'")

        # Read data from all seeds for this basename
        all_data = []
        for seed_dir in seed_dirs:
            result = read_infidelity(seed_dir)
            if result is not None:
                all_data.append(result)

        if len(all_data) == 0:
            print(f"Warning: No valid data found for basename '{basename}'")
            continue

        print(f"Successfully read data from {len(all_data)} seeds for '{basename}'")

        # Filter out traces that stop at infidelity > 1.2e-6
        filtered_data = []
        threshold = 1.2e-6
        for iteration, infidelity, seed_num, x0_norm in all_data:
            final_infidelity = infidelity[-1]
            if final_infidelity <= threshold:
                filtered_data.append((iteration, infidelity, seed_num, x0_norm))
            else:
                print(f"  Excluding seed {seed_num}: stopped at infidelity {final_infidelity:.2e} > {threshold:.2e}")

        if len(filtered_data) == 0:
            print(f"Warning: No seeds reached infidelity <= {threshold:.2e} for basename '{basename}'")
            continue

        print(f"Using {len(filtered_data)} seeds that reached infidelity <= {threshold:.2e} (excluded {len(all_data) - len(filtered_data)})")

        # Use filtered data from here on
        all_data = filtered_data

        # Create label for this dataset
        if len(basenames) > 1:
            dataset_label = basename.replace('_seed', '').replace('config_', '')
        else:
            dataset_label = None

        # Plot individual seeds
        if show_individual:
            for iteration, infidelity, seed_num, x0_norm in all_data:
                ax.plot(iteration, infidelity,
                       linestyle='-', linewidth=0.8, alpha=0.3,
                       color=color, label='_nolegend_')

            # Add one labeled line for the legend
            if dataset_label:
                label = f'{dataset_label} individuals (n={len(all_data)})'
            else:
                label = f'Individual seeds (n={len(all_data)})'
            ax.plot([], [], linestyle='-', linewidth=0.8, alpha=0.3,
                   color=color, label=label)

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

            for iteration, infidelity, seed_num, x0_norm in all_data:
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

                # Plot mean with dataset-specific color
                if dataset_label:
                    mean_label = f'{dataset_label} mean'
                    std_label = f'{dataset_label} ±1σ'
                else:
                    mean_label = 'Mean'
                    std_label = '±1 std dev'

                ax.plot(mean_iter_valid, infidelity_valid,
                       linestyle='-', linewidth=2.5, color=color,
                       marker='o', markersize=5, markevery=max(1, len(mean_iter_valid)//20),
                       label=mean_label, zorder=10)

                # Plot standard deviation band
                ax.fill_betweenx(infidelity_valid,
                               mean_iter_valid - std_iter_valid,
                               mean_iter_valid + std_iter_valid,
                               alpha=0.2, color=color, label=std_label)

                # Also show min and max
                if dataset_label:
                    minmax_label = f'{dataset_label} min/max'
                else:
                    minmax_label = 'Min/Max'

                ax.plot(min_iter_valid, infidelity_valid,
                       linestyle='--', linewidth=1.5, color=color,
                       alpha=0.4, label=minmax_label, zorder=5)
                ax.plot(max_iter_valid, infidelity_valid,
                       linestyle='--', linewidth=1.5, color=color,
                       alpha=0.4, label='_nolegend_', zorder=5)

    # Format plot
    # Get directory name for title
    dir_name = os.path.basename(os.path.abspath(parent_dir))

    # Build title with optional suffix
    if title_suffix:
        title = f'Infidelity Convergence: {dir_name} {title_suffix}'
    else:
        title = f'Infidelity Convergence: {dir_name}'

    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Infidelity', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc='best')
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
    # Create filename based on basenames
    if len(basenames) == 1:
        basename_str = basenames[0].replace('_seed', '')
    else:
        basename_str = 'comparison'

    if title_suffix:
        # Create filename-safe version of suffix (replace spaces and special chars)
        suffix_safe = title_suffix.replace(' ', '_').replace('/', '_').replace('\\', '_')
        output_filename = os.path.join(parent_dir, f'seeds_infidelity_{basename_str}_{suffix_safe}.png')
    else:
        output_filename = os.path.join(parent_dir, f'seeds_infidelity_{basename_str}.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_filename}")

    # Show the plot
    print("\nDisplaying plot... Close the window to exit.")
    plt.show(block=True)

def plot_convergence_vs_initial_norm(parent_dir='.', basenames=None, title_suffix=''):
    """
    Create scatter plot of iterations to convergence vs initial control vector norm.

    Parameters
    ----------
    parent_dir : str
        Parent directory containing seed subdirectories
    basenames : list of str
        List of basename patterns for seed directories (default: ['seed'])
    title_suffix : str
        Optional suffix to add to title and filename (default: '')
    """
    if basenames is None:
        basenames = ['seed']

    # Ensure basenames is a list
    if isinstance(basenames, str):
        basenames = [basenames]

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Color cycle for different basenames
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

    # Process each basename separately
    for basename_idx, basename in enumerate(basenames):
        color = colors[basename_idx % len(colors)]

        # Find all seed directories for this basename
        seed_dirs = find_seed_directories(parent_dir, basename)

        if len(seed_dirs) == 0:
            print(f"Warning: No {basename}_* directories found in {parent_dir}")
            continue

        print(f"Found {len(seed_dirs)} directories for basename '{basename}'")

        # Read data from all seeds for this basename
        all_data = []
        for seed_dir in seed_dirs:
            result = read_infidelity(seed_dir)
            if result is not None:
                all_data.append(result)

        if len(all_data) == 0:
            print(f"Warning: No valid data found for basename '{basename}'")
            continue

        print(f"Successfully read data from {len(all_data)} seeds for '{basename}'")

        # Filter out traces that stop at infidelity > 1.2e-6 and extract convergence data
        threshold = 1.2e-6
        x0_norms = []
        xf_norms = []
        iters_to_converge = []

        for iteration, infidelity, seed_num, x0_norm in all_data:
            final_infidelity = infidelity[-1]
            if final_infidelity <= threshold and x0_norm is not None:
                # Read final control norm from the same file
                seed_dir = seed_dirs[seed_num - 1] if seed_num <= len(seed_dirs) else None
                if seed_dir:
                    history_file = os.path.join(seed_dir, 'optim_history.dat')
                    try:
                        data = np.loadtxt(history_file, skiprows=1)
                        if data.ndim == 1:
                            data = data.reshape(1, -1)
                        xf_norm = data[-1, 12] if data.shape[1] >= 13 else None
                        if xf_norm is not None:
                            x0_norms.append(x0_norm)
                            xf_norms.append(xf_norm)
                            iters_to_converge.append(iteration[-1])  # Last iteration number
                        else:
                            print(f"  Excluding seed {seed_num}: no final norm data")
                    except:
                        print(f"  Excluding seed {seed_num}: error reading final norm")
            else:
                if x0_norm is None:
                    print(f"  Excluding seed {seed_num}: no initial norm data")
                else:
                    print(f"  Excluding seed {seed_num}: stopped at infidelity {final_infidelity:.2e} > {threshold:.2e}")

        if len(x0_norms) == 0:
            print(f"Warning: No seeds with valid convergence data for basename '{basename}'")
            continue

        print(f"Using {len(x0_norms)} seeds with valid convergence data (excluded {len(all_data) - len(x0_norms)})")

        # Create label for this dataset
        if len(basenames) > 1:
            dataset_label = basename.replace('_seed', '').replace('config_', '')
        else:
            dataset_label = None

        # Scatter plot - initial and final norms
        if dataset_label:
            label_init = f'{dataset_label} initial (n={len(x0_norms)})'
            label_final = f'{dataset_label} final'
        else:
            label_init = f'Initial norm (n={len(x0_norms)})'
            label_final = f'Final norm'

        # Plot initial norms with circles
        ax.scatter(x0_norms, iters_to_converge,
                  s=80, alpha=0.6, color=color, label=label_init,
                  marker='o', edgecolors='black', linewidth=0.5)

        # Plot final norms with squares
        ax.scatter(xf_norms, iters_to_converge,
                  s=80, alpha=0.6, color=color, label=label_final,
                  marker='s', edgecolors='black', linewidth=0.5)

    # Format plot
    # Get directory name for title
    dir_name = os.path.basename(os.path.abspath(parent_dir))

    # Build title with optional suffix
    if title_suffix:
        title = f'Convergence vs Initial Control Norm: {dir_name} {title_suffix}'
    else:
        title = f'Convergence vs Initial Control Norm: {dir_name}'

    ax.set_xlabel('Control Vector Norm ||x||/√D', fontsize=12)
    ax.set_ylabel('Iterations to Convergence', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc='best', ncol=2)

    plt.tight_layout()

    # Save the figure
    # Create filename based on basenames
    if len(basenames) == 1:
        basename_str = basenames[0].replace('_seed', '')
    else:
        basename_str = 'comparison'

    if title_suffix:
        # Create filename-safe version of suffix (replace spaces and special chars)
        suffix_safe = title_suffix.replace(' ', '_').replace('/', '_').replace('\\', '_')
        output_filename = os.path.join(parent_dir, f'convergence_vs_x0norm_{basename_str}_{suffix_safe}.png')
    else:
        output_filename = os.path.join(parent_dir, f'convergence_vs_x0norm_{basename_str}.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Scatter plot saved to: {output_filename}")

    # Show the plot
    print("\nDisplaying scatter plot... Close the window to exit.")
    plt.show(block=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plot infidelity convergence from multiple seed directories.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                                          # Default: current directory, seed_*
  %(prog)s qft2                                                     # Specify parent directory
  %(prog)s qft2 --suffix "MINRES+LR"                                # Add title suffix
  %(prog)s . --basename "config_QFT2_dual_seed"                     # Single basename
  %(prog)s . --basename "config_QFT2_lbfgs_seed" --suffix "L-BFGS"  # Single basename with suffix
  %(prog)s . --basename "config_QFT2_dual_seed" "config_QFT2_lbfgs_seed"  # Compare two datasets

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
        nargs='+',
        default=['seed'],
        help='Basename pattern(s) for seed directories (default: "seed", looks for seed_*). Can specify multiple for comparison.'
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

    # Plot infidelity convergence
    plot_seeds_infidelity(
        args.directory,
        args.basename,
        args.suffix,
        show_individual=not args.no_individual,
        show_statistics=not args.no_statistics
    )

    # Plot convergence vs initial norm scatter plot
    plot_convergence_vs_initial_norm(
        args.directory,
        args.basename,
        args.suffix
    )
