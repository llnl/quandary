#!/usr/bin/env python3
"""
Visualize KSP convergence history from multiple ksp_convergence_*.dat files

The script can plot multiple convergence histories on the same figure for comparison.
It ALWAYS saves the plot to a PNG file and optionally displays it interactively.

Usage:
    python3 plot_ksp_convergence.py [file1.dat file2.dat ...] [--output output.png]

Arguments:
    file1.dat file2.dat  : One or more KSP convergence data files (optional)
    --output             : Output plot filename (optional)

If no input files are provided, the script will auto-detect all ksp_convergence_*.dat
files in the current directory and plot them together.

Plots residual norm and solution norm as functions of iteration number
using log scale for the y-axis. Each file gets a different color/marker.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
import argparse

def read_ksp_file(filename):
    """
    Read KSP convergence data and metadata from file.

    Returns
    -------
    dict with keys: ksp_type, converged_reason, rtol, maxiter, iterations, residual_norm, solution_norm
    """
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found.")
        return None

    # Read metadata from comment lines
    ksp_type = "Unknown"
    converged_reason = "Unknown"
    rtol = None
    maxiter = None
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('# KSP_TYPE:'):
                ksp_type = line.split(':', 1)[1].strip()
            elif line.startswith('# CONVERGED_REASON:'):
                converged_reason = line.split(':', 1)[1].strip()
            elif line.startswith('# RTOL:'):
                rtol = float(line.split(':', 1)[1].strip())
            elif line.startswith('# MAXITER:'):
                maxiter = int(line.split(':', 1)[1].strip())

    # Read data, skipping comment lines
    data = np.loadtxt(filename, comments='#')

    if data.size == 0:
        print(f"Warning: No data found in '{filename}'.")
        return None

    # Extract columns
    iterations = data[:, 0].astype(int)
    residual_norm = data[:, 1]
    solution_norm = data[:, 2]

    return {
        'filename': filename,
        'ksp_type': ksp_type,
        'converged_reason': converged_reason,
        'rtol': rtol,
        'maxiter': maxiter,
        'iterations': iterations,
        'residual_norm': residual_norm,
        'solution_norm': solution_norm
    }


def plot_ksp_convergence_multi(filenames, output_file=None, show_interactive=False):
    """
    Plot KSP convergence history from multiple data files.

    Parameters
    ----------
    filenames : list of str
        Paths to KSP convergence data files
    output_file : str
        Path for the output plot image
    show_interactive : bool
        Whether to display the plot interactively
    """
    # Read all files
    datasets = []
    for filename in filenames:
        data = read_ksp_file(filename)
        if data is not None:
            datasets.append(data)

    if len(datasets) == 0:
        print("Error: No valid data files to plot.")
        sys.exit(1)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Color and marker cycling
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    linestyles = ['-', '--', '-.', ':']

    # Plot each dataset
    for idx, data in enumerate(datasets):
        color = colors[idx % len(colors)]
        marker = markers[idx % len(markers)]
        linestyle = linestyles[idx % len(linestyles)]

        # Extract iteration number from filename if present (format: ksp_convergence_TYPE_iterNNNN.dat)
        import re
        filename = os.path.basename(data['filename'])
        iter_match = re.search(r'_iter(\d+)', filename)

        ksp_label = data['ksp_type'].upper()
        if iter_match:
            iter_num = int(iter_match.group(1))
            label = f"Iter {iter_num} - {ksp_label} (n={len(data['iterations'])})"
        else:
            label = f"{ksp_label} (n={len(data['iterations'])})"

        # Plot residual norm
        ax1.semilogy(data['iterations'], data['residual_norm'],
                     marker=marker, linestyle=linestyle, linewidth=2, markersize=6,
                     color=color, label=label, markevery=max(1, len(data['iterations'])//10))

        # Plot solution norm
        ax2.semilogy(data['iterations'], data['solution_norm'],
                     marker=marker, linestyle=linestyle, linewidth=2, markersize=6,
                     color=color, label=label, markevery=max(1, len(data['iterations'])//10))

    # Configure residual norm plot
    ax1.set_xlabel('Iteration', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Residual Norm', fontsize=12, fontweight='bold')
    ax1.set_title('KSP Residual Convergence', fontsize=14, fontweight='bold')
    ax1.grid(True, which='both', alpha=0.3, linestyle='--')
    ax1.legend(fontsize=10, loc='best')

    # Configure solution norm plot
    ax2.set_xlabel('Iteration', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Solution Norm', fontsize=12, fontweight='bold')
    ax2.set_title('KSP Solution Norm Evolution', fontsize=14, fontweight='bold')
    ax2.grid(True, which='both', alpha=0.3, linestyle='--')
    ax2.legend(fontsize=10, loc='best')

    # Adjust layout
    plt.tight_layout()

    # Generate output filename if not provided
    if output_file is None:
        if len(datasets) == 1:
            base_name = os.path.splitext(datasets[0]['filename'])[0]
            output_file = f"{base_name}.png"
        else:
            output_file = "ksp_convergence_comparison.png"

    # Always save to file
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")

    # Show interactively if requested
    if show_interactive:
        print("Displaying plot interactively (close window to continue)...")
        plt.show(block=True)

    # Print summary statistics for each dataset
    print(f"\n{'='*70}")
    print(f"KSP Convergence Summary ({len(datasets)} dataset(s)):")
    print(f"{'='*70}")
    for idx, data in enumerate(datasets):
        print(f"\n[{idx+1}] {os.path.basename(data['filename'])}")
        print(f"  KSP type: {data['ksp_type']}")
        print(f"  Convergence reason: {data['converged_reason']}")
        if data['rtol'] is not None:
            print(f"  Relative tolerance (rtol): {data['rtol']:.6e}")
        if data['maxiter'] is not None:
            print(f"  Maximum iterations: {data['maxiter']}")
        print(f"  Total iterations: {len(data['iterations'])}")
        print(f"  Initial residual norm: {data['residual_norm'][0]:.6e}")
        print(f"  Final residual norm: {data['residual_norm'][-1]:.6e}")
        print(f"  Reduction factor: {data['residual_norm'][0]/data['residual_norm'][-1]:.2e}")
        print(f"  Final solution norm: {data['solution_norm'][-1]:.6e}")
    print(f"{'='*70}")


if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Visualize KSP convergence history from multiple files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 plot_ksp_convergence.py                                        # Auto-detect all *.dat files
  python3 plot_ksp_convergence.py ksp_convergence_cg.dat                # Plot single file
  python3 plot_ksp_convergence.py ksp_*.dat                             # Plot all matching files
  python3 plot_ksp_convergence.py cg.dat minres.dat --output comp.png   # Compare two files
        """
    )
    parser.add_argument('input_files', nargs='*', default=[],
                        help='KSP convergence data files (e.g., ksp_convergence_cg.dat ksp_convergence_minres.dat)')
    parser.add_argument('--output', '-o', dest='output_file', default=None,
                        help='Output plot filename (if omitted, auto-generates and displays interactively)')

    args = parser.parse_args()

    # Determine input files
    input_files = args.input_files
    if not input_files:
        # Try to auto-detect ksp_convergence_*.dat files
        pattern = 'ksp_convergence_*.dat'
        input_files = sorted(glob.glob(pattern))
        if input_files:
            print(f"Auto-detected {len(input_files)} file(s): {', '.join(input_files)}")
        else:
            # Fall back to old filename for backward compatibility
            if os.path.exists('ksp_convergence.dat'):
                input_files = ['ksp_convergence.dat']
                print(f"Using: ksp_convergence.dat")
            else:
                print("Error: No KSP convergence data file found.")
                print("Please provide input files or run the simulation first.")
                sys.exit(1)

    # Check all files exist
    for f in input_files:
        if not os.path.exists(f):
            print(f"Error: File '{f}' not found.")
            sys.exit(1)

    # Determine output file and whether to show interactively
    output_file = args.output_file
    show_interactive = (args.output_file is None)  # Show interactively if no output file specified

    # Plot the convergence history
    plot_ksp_convergence_multi(input_files, output_file, show_interactive)
