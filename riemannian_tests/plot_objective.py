#!/usr/bin/env python3
"""
Plot the Objective function convergence from optimization history file.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_objective(filenames):
    """
    Plot objective function vs iteration number for one or more files.

    Parameters
    ----------
    filenames : str or list of str
        Path(s) to the data file(s) (format: optim_history_*.dat)
    """
    # Ensure filenames is a list
    if isinstance(filenames, str):
        filenames = [filenames]

    # Create figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 14))

    # First pass: find the maximum number of data points across all files
    max_points = 0
    for filename in filenames:
        data = np.loadtxt(filename, skiprows=1)
        max_points = max(max_points, len(data))

    # Calculate marker frequency based on the longest dataset (show ~20 markers)
    marker_every = max(1, max_points // 20)

    # Plot each file
    for idx, filename in enumerate(filenames):
        # Read the data, skipping the header line
        data = np.loadtxt(filename, skiprows=1)

        # Extract columns
        iteration = data[:, 0].astype(int)
        objective = data[:, 1]
        infidelity = data[:, 5]  # Column 6 (1-indexed) = index 5 (0-indexed)
        tikhonov = data[:, 6]    # Column 7 (1-indexed) = index 6 (0-indexed)
        ksp_iters = data[:, 11].astype(int)

        # Extract label from filename (remove path and extension)
        import os
        label = os.path.basename(filename).replace('.dat', '').replace('optim_history_', '')

        # Top subplot: Objective and Tikhonov (log scale)
        ax1.plot(iteration, objective, 'o-', linewidth=2, markersize=4,
                label=f'{label} (Obj)', color=f'C{idx}', markevery=marker_every)
        ax1.plot(iteration, tikhonov, 's--', linewidth=1.5, markersize=3,
                label=f'{label} (Tikh)', color=f'C{idx}', alpha=0.7, markevery=marker_every)

        # Middle subplot: Infidelity (log scale)
        # Use alternating line styles and semi-transparency for better distinguishability
        linestyles = ['-', '--', '-.', ':']
        ax2.plot(iteration, infidelity, marker='o', linestyle=linestyles[idx % len(linestyles)],
                linewidth=2, markersize=4, alpha=0.7,
                label=label, color=f'C{idx}', markevery=marker_every)

        # Bottom subplot: KSP iterations (linear scale)
        # Use alternating line styles and semi-transparency for better distinguishability
        ax3.plot(iteration, ksp_iters, marker='o', linestyle=linestyles[idx % len(linestyles)],
                linewidth=2, markersize=4, alpha=0.7,
                label=label, color=f'C{idx}', markevery=marker_every)

    # Format top subplot
    ax1.set_xlabel('Iteration', fontsize=12)
    ax1.set_ylabel('Objective', fontsize=12)
    ax1.set_title('Objective Function Convergence', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_yscale('log')

    # Format middle subplot
    ax2.set_xlabel('Iteration', fontsize=12)
    ax2.set_ylabel('Infidelity', fontsize=12)
    ax2.set_title('Infidelity Convergence', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_yscale('log')

    # Format bottom subplot
    ax3.set_xlabel('Iteration', fontsize=12)
    ax3.set_ylabel('KSP iterations', fontsize=12)
    ax3.set_title('KSP Iterations per Gauss-Newton Step', fontsize=14)
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=10)

    plt.tight_layout()

    # Save the figure
    if len(filenames) == 1:
        output_filename = filenames[0].replace('.dat', '_objective.png')
    else:
        output_filename = 'comparison_objective.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_filename}")

    # Show the plot interactively and wait for window to close
    print("Displaying plot... Close the window to exit.")
    plt.show(block=True)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        # Default filename
        filenames = ["optim_history_ksp_damping.dat"]
    else:
        # Accept multiple filenames
        filenames = sys.argv[1:]

    plot_objective(filenames)
