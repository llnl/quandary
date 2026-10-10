#!/usr/bin/env python3
"""
Plot Gauss-Newton matrix as a heatmap in log scale.

Usage:
    python plot_gn_matrix.py <matrix_file.dat>

The script reads a PETSc binary matrix file and displays it as a log-scale heatmap.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def read_petsc_matrix(filename):
    """Read a PETSc binary matrix file using PETSc's built-in Python module."""

    # Try to find PetscBinaryIO module in PETSc installation
    petsc_dir = os.environ.get('PETSC_DIR')
    if not petsc_dir:
        # Try common homebrew location
        petsc_dir = '/opt/homebrew/Cellar/petsc/3.25.5'

    petsc_python_path = os.path.join(petsc_dir, 'lib', 'petsc', 'bin')

    if not os.path.exists(petsc_python_path):
        # Try alternative path
        petsc_python_path = os.path.join(petsc_dir, 'share', 'petsc', 'lib', 'petsc', 'bin')

    if os.path.exists(petsc_python_path):
        sys.path.insert(0, petsc_python_path)

    try:
        import PetscBinaryIO

        print(f"Using PetscBinaryIO from: {petsc_python_path}")

        # Read the binary file
        io = PetscBinaryIO.PetscBinaryIO()

        # readBinaryFile returns a list of objects
        objects = io.readBinaryFile(filename)

        if len(objects) == 0:
            print("Error: No objects found in file")
            sys.exit(1)

        # Get the first object (the dense matrix saved by Quandary)
        mat = objects[0]

        # MatSparse is a 2-element tuple: ((nrows, ncols), (row_offsets, col_indices, values))
        # This is CSR (Compressed Sparse Row) format
        if hasattr(mat, '__class__') and 'MatSparse' in mat.__class__.__name__:
            # Extract dimensions and CSR data
            (nrows, ncols), (row_offsets, col_indices, values) = mat

            print(f"Matrix dimensions: {nrows} x {ncols}")
            print(f"Number of stored entries: {len(values)}")

            # Convert CSR format to dense
            dense = np.zeros((nrows, ncols))

            # For each row
            for i in range(nrows):
                # Get the range of indices for this row
                start = row_offsets[i]
                end = row_offsets[i + 1]

                # Fill in the values for this row
                for idx in range(start, end):
                    j = col_indices[idx]
                    dense[i, j] = values[idx]

            print(f"Successfully converted CSR format to dense")
            return dense

        elif isinstance(mat, np.ndarray):
            # Already a numpy array
            print(f"Successfully read matrix: {mat.shape[0]} x {mat.shape[1]}")
            return mat
        else:
            print(f"Error: Unexpected matrix format: {type(mat)}")
            sys.exit(1)

    except ImportError as e:
        print(f"Error: Could not import PetscBinaryIO: {e}")
        print(f"Searched in: {petsc_python_path}")
        print("\nAlternative: Set PETSC_DIR environment variable to your PETSc installation")
        print("Example: export PETSC_DIR=/opt/homebrew/Cellar/petsc/3.25.5")
        sys.exit(1)


def plot_matrix_heatmap(matrix, filename, output_file=None):
    """Plot matrix as a log-scale heatmap."""

    # Get matrix properties
    n, m = matrix.shape
    print(f"Matrix dimensions: {n} x {m}")
    print(f"Matrix min value: {np.min(matrix):.6e}")
    print(f"Matrix max value: {np.max(matrix):.6e}")
    print(f"Matrix mean value: {np.mean(matrix):.6e}")

    # Check if matrix is all zeros
    max_val = np.max(np.abs(matrix))
    if max_val == 0:
        print("\n*** WARNING: Matrix is all zeros! ***")
        print("The Gauss-Newton matrix was not properly computed.")
        print("Check that xeval_GN is set correctly before matrix assembly.\n")
        return

    # Check if matrix is symmetric
    if n == m:
        mat_norm = np.linalg.norm(matrix)
        if mat_norm > 0:
            symmetry_error = np.linalg.norm(matrix - matrix.T) / mat_norm
            print(f"Symmetry error: {symmetry_error:.6e}")
            if symmetry_error < 1e-10:
                print("Matrix is symmetric (within numerical tolerance)")
        else:
            print("Cannot check symmetry (zero matrix)")

    # Take absolute values for log scale (handles negative values)
    mat_abs = np.abs(matrix)

    # Set very small values to a minimum threshold for visualization
    min_threshold = np.max(mat_abs) * 1e-16
    mat_abs[mat_abs < min_threshold] = min_threshold

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot with log scale (origin='upper' puts row 0 at top, standard matrix convention)
    im = ax.imshow(mat_abs, cmap='viridis', norm=LogNorm(), aspect='auto', origin='upper')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('|Matrix elements| (log scale)', rotation=270, labelpad=20)

    # Labels and title
    ax.set_xlabel('Column index')
    ax.set_ylabel('Row index')

    # Determine matrix type from filename
    import os
    base_filename = os.path.basename(filename)

    if 'dual' in filename.lower():
        title = f'Dual Gauss-Newton Matrix (L L^T)\n{base_filename}\nDimensions: {n} x {m}'
    elif 'primal' in filename.lower():
        title = f'Primal Gauss-Newton Matrix (L^T L)\n{base_filename}\nDimensions: {n} x {m}'
    else:
        title = f'Gauss-Newton Matrix\n{base_filename}\nDimensions: {n} x {m}'

    ax.set_title(title)

    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()

    # Save to file
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_file}")

    # Also show interactively
    plt.show(block=True)
    plt.close()


def plot_eigenvalue_distribution(matrix, output_file=None, num_top_eigs=20, show_top_eigs=False, input_filename=""):
    """Plot eigenvalue distribution of the matrix."""

    if matrix.shape[0] != matrix.shape[1]:
        print("Warning: Matrix is not square, skipping eigenvalue plot")
        return

    # Check if matrix is all zeros
    if np.max(np.abs(matrix)) == 0:
        print("Skipping eigenvalue computation (zero matrix)")
        return

    print("\nComputing eigenvalues...")
    eigenvalues = np.linalg.eigvalsh(matrix)  # For symmetric matrices
    eigenvalues = np.sort(eigenvalues)[::-1]  # Sort descending

    print(f"Largest eigenvalue: {eigenvalues[0]:.6e}")
    print(f"Smallest eigenvalue: {eigenvalues[-1]:.6e}")
    print(f"Condition number: {eigenvalues[0] / eigenvalues[-1]:.6e}")

    # Get base filename for title
    import os
    base_filename = os.path.basename(input_filename)

    # Ensure we don't try to plot more eigenvalues than exist
    num_top_eigs_clamped = min(num_top_eigs, len(eigenvalues))

    # Show separate subplot for top eigenvalues if user explicitly requested it
    if show_top_eigs:
        # For primal matrices, create two subplots: full spectrum and top 20
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Plot 1: Full eigenvalue spectrum
        ax1.semilogy(range(len(eigenvalues)), eigenvalues, 'b.-', linewidth=1, markersize=3)
        ax1.set_xlabel('Eigenvalue index')
        ax1.set_ylabel('Eigenvalue (log scale)')
        ax1.set_title(f'Full Eigenvalue Spectrum\n{base_filename}')
        ax1.grid(True, alpha=0.3)

        # Plot 2: First N eigenvalues (clamped to available)
        ax2.semilogy(range(num_top_eigs_clamped), eigenvalues[:num_top_eigs_clamped], 'r.-', linewidth=1, markersize=3)
        ax2.set_xlabel('Eigenvalue index')
        ax2.set_ylabel('Eigenvalue (log scale)')
        if num_top_eigs_clamped < num_top_eigs:
            ax2.set_title(f'First {num_top_eigs_clamped} Eigenvalues (all available)\n{base_filename}')
        else:
            ax2.set_title(f'First {num_top_eigs_clamped} Eigenvalues\n{base_filename}')
        ax2.grid(True, alpha=0.3)
    else:
        # For dual or small matrices, show only full spectrum
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot eigenvalue spectrum
        ax.semilogy(range(len(eigenvalues)), eigenvalues, 'b.-', linewidth=1, markersize=3)
        ax.set_xlabel('Eigenvalue index')
        ax.set_ylabel('Eigenvalue (log scale)')
        ax.set_title(f'Eigenvalue Spectrum\n{base_filename}')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save to file (derive eigenvalue filename from heatmap filename)
    base = output_file.rsplit('.', 1)[0].replace('_heatmap', '')
    eig_file = base + '_eigenvalues.png'
    plt.savefig(eig_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {eig_file}")

    # Also show interactively
    plt.show(block=True)
    plt.close()


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_gn_matrix.py <matrix_file.dat> [output_file_prefix.png] [num_top_eigenvalues]")
        print("\nGenerates two plots:")
        print("  1. Matrix heatmap (log scale)")
        print("  2. Eigenvalue distribution")
        print("\nFor primal matrices, the eigenvalue plot shows both full spectrum and top N eigenvalues.")
        print("\nExamples:")
        print("  python plot_gn_matrix.py GN_dual_matrix.dat")
        print("    -> Creates GN_dual_matrix_heatmap.png and GN_dual_matrix_eigenvalues.png")
        print("\n  python plot_gn_matrix.py GN_primal_matrix.dat my_plot.png 30")
        print("    -> Creates my_plot.png and my_plot_eigenvalues.png, showing top 30 eigenvalues")
        sys.exit(1)

    input_file = sys.argv[1]

    # Smart argument parsing: if argv[2] is a number, treat it as num_top_eigs
    num_top_eigs = 20
    output_file = None
    show_top_eigs = False  # Track if num_top_eigs was explicitly provided

    if len(sys.argv) > 2:
        # Try to parse second argument as integer
        try:
            num_top_eigs = int(sys.argv[2])
            if num_top_eigs <= 0:
                print("Error: num_top_eigenvalues must be positive")
                sys.exit(1)
            # It's a number, so use default output filename
            output_file = None
            show_top_eigs = True  # User explicitly provided number
        except ValueError:
            # Not a number, treat as output filename
            output_file = sys.argv[2]

    if len(sys.argv) > 3:
        # Third argument is num_top_eigs (if second was filename)
        try:
            num_top_eigs = int(sys.argv[3])
            if num_top_eigs <= 0:
                print("Error: num_top_eigenvalues must be positive")
                sys.exit(1)
            show_top_eigs = True  # User explicitly provided number
        except ValueError:
            print("Error: Third argument (num_top_eigenvalues) must be an integer")
            sys.exit(1)

    # Generate default output filename if not provided
    if output_file is None:
        base = input_file.rsplit('.', 1)[0]
        output_file = base + '_heatmap.png'

    print(f"Reading matrix from: {input_file}")
    matrix = read_petsc_matrix(input_file)

    print("\n=== Matrix Heatmap ===")
    plot_matrix_heatmap(matrix, input_file, output_file)

    print("\n=== Eigenvalue Analysis ===")
    plot_eigenvalue_distribution(matrix, output_file, num_top_eigs, show_top_eigs, input_file)


if __name__ == "__main__":
    main()
