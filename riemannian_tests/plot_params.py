#!/usr/bin/env python
"""
Plot parameter data as scatter plots showing all data points with min/max values indicated.

The script reads params.dat files and creates scatter plots of the parameter values.

Usage:
    python3 plot_params.py [file1.dat file2.dat ...] [--output output.png]

Arguments:
    file1.dat file2.dat  : One or more parameter data files (optional)
    --output             : Output plot filename (optional)

If no input files are provided, the script will auto-detect all params*.dat
files in the current directory and plot them together.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
import argparse

def plot_params(filenames, output_file=None):
    """
    Plot parameter data from one or more params data files.

    Parameters
    ----------
    filenames : list of str
        Paths to parameter data files
    output_file : str
        Path for the output plot image (optional)
    """
    if isinstance(filenames, str):
        filenames = [filenames]

    # Create figure with subplots (one per file)
    n_files = len(filenames)
    fig, axes = plt.subplots(n_files, 1, figsize=(12, 6*n_files))

    # Make axes iterable even if there's only one file
    if n_files == 1:
        axes = [axes]

    # Plot each file in its own subplot
    for idx, (filename, ax) in enumerate(zip(filenames, axes)):
        # Read the data
        try:
            # Try reading with header first
            data = np.loadtxt(filename, skiprows=1)
        except:
            # If that fails, read without skipping
            data = np.loadtxt(filename)

        # Handle different data shapes
        if data.ndim == 1:
            # Single column - use index as x-axis
            indices = np.arange(len(data))
            values = data
        else:
            # Multiple columns - first column is index/iteration
            indices = data[:, 0]
            values = data[:, 1] if data.shape[1] > 1 else data[:, 0]

        # Scale values by 1/(2*pi) to convert to GHz
        values = values / (2 * np.pi)

        # Calculate min and max
        min_val = np.min(values)
        max_val = np.max(values)
        min_idx = indices[np.argmin(values)]
        max_idx = indices[np.argmax(values)]

        # Extract label from filename (e.g., params.dat -> Parameters)
        basename = os.path.basename(filename)
        label_base = basename.replace('.dat', '').replace('params', 'Parameters')

        # Create scatter plot
        ax.scatter(indices, values, s=50, alpha=0.6, color='C0',
                  edgecolors='black', linewidth=0.5, label='Parameter values')

        # Mark min and max points
        ax.scatter([min_idx], [min_val], s=150, marker='v',
                  color='red', edgecolors='darkred', linewidth=2,
                  label=f'Min = {min_val:.6e} @ {min_idx:.0f}', zorder=10)
        ax.scatter([max_idx], [max_val], s=150, marker='^',
                  color='green', edgecolors='darkgreen', linewidth=2,
                  label=f'Max = {max_val:.6e} @ {max_idx:.0f}', zorder=10)

        # Add horizontal lines at min/max values
        ax.axhline(y=min_val, color='red', linestyle='--',
                  linewidth=1, alpha=0.5)
        ax.axhline(y=max_val, color='green', linestyle='--',
                  linewidth=1, alpha=0.5)

        # Format subplot
        ax.set_xlabel('Index / Iteration', fontsize=12)
        ax.set_ylabel('Amplitude [GHz]', fontsize=12)
        ax.set_title(f'{label_base}', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc='best')

        # Print statistics to console
        mean_val = np.mean(values)
        std_val = np.std(values)
        print(f"\n{label_base} statistics:")
        print(f"  Min:  {min_val:.6e} at index {min_idx:.0f}")
        print(f"  Max:  {max_val:.6e} at index {max_idx:.0f}")
        print(f"  Mean: {mean_val:.6e}")
        print(f"  Std:  {std_val:.6e}")
        print(f"  Range: {max_val - min_val:.6e}")

    plt.tight_layout()

    # Save the figure
    if output_file is None:
        if len(filenames) == 1:
            output_file = filenames[0].replace('.dat', '_plot.png')
        else:
            output_file = 'params_comparison.png'

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")

    # Show the plot interactively and wait for window to close
    print("Displaying plot... Close the window to exit.")
    plt.show(block=True)


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Plot parameter values as scatter plot with min/max indicators',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 plot_params.py                                # Auto-detect all params*.dat files
  python3 plot_params.py params.dat                    # Plot single file
  python3 plot_params.py params*.dat                   # Plot all matching files
  python3 plot_params.py params0.dat params1.dat       # Compare multiple files
  python3 plot_params.py params.dat --output plot.png  # Specify output filename
        """
    )
    parser.add_argument('input_files', nargs='*', default=[],
                        help='Parameter data files (e.g., params.dat params0.dat)')
    parser.add_argument('--output', '-o', dest='output_file', default=None,
                        help='Output plot filename (if omitted, auto-generates)')

    args = parser.parse_args()

    # Determine input files
    input_files = args.input_files
    if not input_files:
        # Try to auto-detect params*.dat files
        pattern = 'params*.dat'
        input_files = sorted(glob.glob(pattern))
        if input_files:
            print(f"Auto-detected {len(input_files)} file(s): {', '.join(input_files)}")
        else:
            print("Error: No params data file found.")
            print("Please provide input files or run the simulation first.")
            sys.exit(1)

    # Check all files exist
    for f in input_files:
        if not os.path.exists(f):
            print(f"Error: File '{f}' not found.")
            sys.exit(1)

    # Plot the parameter data
    plot_params(input_files, args.output_file)
