#!/usr/bin/env python3
"""
Reformat convergence_raw.dat to ksp_convergence format.

Converts from:
  KSP it 0: residual norm = 1.234e-05, solution norm = 0.0

To:
  # KSP convergence history
  # KSP_TYPE: formatted
  # iteration  residual_norm  solution_norm
  0 1.234e-05 0.0

Usage:
  python3 reformat_convergence.py input.dat output.dat
  python3 reformat_convergence.py convergence_raw.dat  # outputs to ksp_convergence_formatted.dat
"""

import sys
import re
import os

def reformat_convergence_file(input_file, output_file=None, ksp_type="formatted"):
    """
    Reformat convergence data from text format to columnar format.

    Parameters
    ----------
    input_file : str
        Input file with "KSP it N: residual norm = X, solution norm = Y" format
    output_file : str, optional
        Output file. If None, generates name from input file.
    ksp_type : str
        KSP solver type to write in header
    """
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)

    # Generate output filename if not provided
    if output_file is None:
        base = os.path.splitext(input_file)[0]
        output_file = f"{base}_formatted.dat"
        if base == "convergence_raw":
            output_file = "ksp_convergence_formatted.dat"

    # Read and parse input file
    data = []
    pattern = r'KSP it (\d+): residual norm = ([0-9.eE+-]+), solution norm = ([0-9.eE+-]+)'

    with open(input_file, 'r') as f:
        for line in f:
            match = re.match(pattern, line)
            if match:
                iteration = int(match.group(1))
                residual = match.group(2)
                solution = match.group(3)
                data.append((iteration, residual, solution))

    if not data:
        print(f"Error: No convergence data found in '{input_file}'")
        sys.exit(1)

    # Write output file
    with open(output_file, 'w') as f:
        f.write("# KSP convergence history\n")
        f.write(f"# KSP_TYPE: {ksp_type}\n")
        f.write("# iteration  residual_norm  solution_norm\n")
        for iteration, residual, solution in data:
            f.write(f"{iteration} {residual} {solution}\n")

    print(f"Reformatted {len(data)} iterations from '{input_file}' to '{output_file}'")
    return output_file


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Reformat convergence data to ksp_convergence format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 reformat_convergence.py convergence_raw.dat
  python3 reformat_convergence.py input.dat output.dat
  python3 reformat_convergence.py raw.dat formatted.dat --ksp-type cg
        """
    )
    parser.add_argument('input_file', help='Input convergence file')
    parser.add_argument('output_file', nargs='?', default=None,
                        help='Output file (default: auto-generated)')
    parser.add_argument('--ksp-type', default='formatted',
                        help='KSP solver type label (default: "formatted")')

    args = parser.parse_args()

    reformat_convergence_file(args.input_file, args.output_file, args.ksp_type)
