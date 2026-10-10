#!/usr/bin/env python3
"""
Plot control parameters (p and q components) as a function of time.

The script reads control*.dat files and plots the drive amplitudes.

Usage:
    python3 plot_control.py [file1.dat file2.dat ...] [--output output.png]

Arguments:
    file1.dat file2.dat  : One or more control data files (optional)
    --output             : Output plot filename (optional)

If no input files are provided, the script will auto-detect all control*.dat
files in the current directory and plot them together.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
import argparse

def plot_control(filenames, output_file=None):
    """
    Plot control parameters from one or more control data files.

    Parameters
    ----------
    filenames : list of str
        Paths to control data files
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
        # Read the data, skipping the header line
        data = np.loadtxt(filename, skiprows=1)

        # Extract columns
        time = data[:, 0]
        drive_p = data[:, 1]
        drive_q = data[:, 2]

        # Extract label from filename (e.g., control0.dat -> Control 0)
        basename = os.path.basename(filename)
        label_base = basename.replace('.dat', '').replace('control', 'Control ')

        # Calculate marker frequency (show ~30 markers regardless of data size)
        n_points = len(time)
        marker_every = max(1, n_points // 30)

        # Plot both p and q components with different colors
        ax.plot(time, drive_p, 'o-', linewidth=2, markersize=4,
                label='p component', color='C0', markevery=marker_every)
        ax.plot(time, drive_q, 's-', linewidth=2, markersize=4,
                label='q component', color='C1', markevery=marker_every)

        # Format subplot
        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Drive amplitude', fontsize=12)
        ax.set_title(f'{label_base}', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc='best')

    plt.tight_layout()

    # Save the figure
    if output_file is None:
        if len(filenames) == 1:
            output_file = filenames[0].replace('.dat', '_plot.png')
        else:
            output_file = 'control_comparison.png'

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")

    # Show the plot interactively and wait for window to close
    print("Displaying plot... Close the window to exit.")
    plt.show(block=True)


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Plot control parameters (p and q) vs time',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 plot_control.py                                   # Auto-detect all control*.dat files
  python3 plot_control.py control0.dat                     # Plot single file
  python3 plot_control.py control*.dat                     # Plot all matching files
  python3 plot_control.py control0.dat control1.dat        # Compare multiple files
  python3 plot_control.py control0.dat --output plot.png  # Specify output filename
        """
    )
    parser.add_argument('input_files', nargs='*', default=[],
                        help='Control data files (e.g., control0.dat control1.dat)')
    parser.add_argument('--output', '-o', dest='output_file', default=None,
                        help='Output plot filename (if omitted, auto-generates)')

    args = parser.parse_args()

    # Determine input files
    input_files = args.input_files
    if not input_files:
        # Try to auto-detect control*.dat files
        pattern = 'control*.dat'
        input_files = sorted(glob.glob(pattern))
        if input_files:
            print(f"Auto-detected {len(input_files)} file(s): {', '.join(input_files)}")
        else:
            print("Error: No control data file found.")
            print("Please provide input files or run the simulation first.")
            sys.exit(1)

    # Check all files exist
    for f in input_files:
        if not os.path.exists(f):
            print(f"Error: File '{f}' not found.")
            sys.exit(1)

    # Plot the control data
    plot_control(input_files, args.output_file)
