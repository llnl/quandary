# Getting Started

## Installation

You have two options to install the Quandary binary:

### Option 1: Install via Spack (Recommended)

If you have [Spack](https://spack.readthedocs.io/en/latest/getting_started.html#installation) installed:
```console
spack install quandary
```

This automatically handles all dependencies including PETSc and MPI.

### Option 2: Build from Source

**Prerequisites for building from source:**

- **PETSc** (Portable, Extensible Toolkit for Scientific Computation)
- **MPI** implementation for parallel execution
- **CMake** 3.23 or later
- **C++ compiler** with C++17 support

**Steps:**

- Clone the repository:

```console
git clone https://github.com/LLNL/quandary.git
cd quandary
```

- Follow the detailed build instructions in the [README.md](https://github.com/LLNL/quandary/blob/main/README.md)

## Running

### C++ Interface

Test your installation with the provided template:
```console
./quandary config_template.toml  # Serial execution
mpirun -np 4 ./quandary config_template.toml  # Parallel execution
```

Results are written as column-based text files in the output directory `data_out/`. The `config_template.toml` is currently set to run a CNOT optimization test case. It lists all available options and configurations, and is filled with comments that should help users to set up new simulation and optimization runs, and match the input options to the equations found in this document.

You can silence Quandary by adding the `--quiet` command line argument.

### Python Interface

For Python interface, after cloning quandary:
```console
pip install -e .
```

Test the Python interface with a working example:
```console
cd examples
python example_cnot.py
```

This example demonstrates:

- Setting up a 2-qubit system
- Defining a CNOT gate target
- Running optimization
- Plotting results

The output will show optimization progress and generate control pulse plots.

## Next Steps

- **[Examples](https://github.com/LLNL/quandary-examples/)**: Check out the quandary-examples repo
- **[Jupyter Example](QuandaryWithPython_HowTo.ipynb)**: Check out the Jupyter Notebook Tutorial in these docs or in the above repo
- **[User Guide](user_guide.md)**: Comprehensive documentation of Quandary's capabilities and details
- **[C++ Config Reference](config.md)**: Reference listing configuration files in detail
- **[Python API Reference](python_api.md)**: Reference listing the Python interface in detail
