# <img src="/quandary_logo/quandary-logo_logo-inline-color.png" width="512" alt="Quandary"/>
[![Build and Test](https://github.com/LLNL/quandary/actions/workflows/test.yml/badge.svg)](https://github.com/LLNL/quandary/actions/workflows/test.yml)

# Optimal control for open and closed quantum systems
Quandary simulates and optimizes the time evolution of closed and open quantum systems, given a Hamiltonian that models driven superconducting quantum devices. The
underlying dynamics are modelled by either Schroedinger's equation (closed systems, state vector), or Lindblad's master equation (open systems, density matrix). Quandary solves the respective ordinary differential equation (ODE) numerically by applying a time-stepping integration scheme, and applies a gradient-based optimization scheme to design optimal control pulses that drive the quantum system to desired targets.
The target can be a unitary gate, i.e. optimizing for pulses that
realize a logical quantum operation, or state preparation that aims to drive the quantum system from one (or multiple) initial state to a desired target state, such as the ground state.

Quandary targets deployment on High-Performance Computing platforms, offering various levels for parallelization using the message passing paradigm. Quandary is written in C++ and executed by providing a text-based configuration file. Further, a Python interface is available to call the underlying C++ code from within a python interpreter, to ease usage.

Feel free to reach out to Stefanie Guenther [guenther5@llnl.gov] for any question you may have. 

# Documentation
Both user and code documentation is available [here](https://software.llnl.gov/quandary/).

# Building
Quandary uses CMake and [BLT](https://github.com/LLNL/blt) to handle builds. Since BLT is included as a
submodule, first make sure you run:
```
git submodule init && git submodule update
```

This project relies on Petsc [https://petsc.org/release/] to handle (parallel) linear algebra. You can either use Spack to install Quandary alongside Petsc, or use CMake to install Quandary given an existing Petsc installation. 

## Building using Spack
Spack can be used to install Quandary, including the required dependency on Petsc, as well as Python packages and interface.

1. To install Spack, clone the repo and add to your shell following the steps [here](https://spack.readthedocs.io/en/latest/getting_started.html#installation).

2. To setup your compilers in your local Spack configuration:
   ```
   spack compiler find
   ```

3. To activate Quandary's spack environment, run the following in your local Quandary folder:
    ```
    spack env activate .spack_env/
    ```

4. Finally, to install the necessary dependencies and build Quandary run:
    ```
    spack install
    ```
    Note: This step could take a while the first time. The second time you run this is should be much faster, only looking for changes in the environment or local code.

Note that `spack install` will build Quandary using CMake from your local source code and install the binary in your Spack environment. To install Quandary's python interface and the python dependencies, use
```
pip install -e .
```

#### Optional Spack environment variations
The Spack environment used to build Quandary is defined in `.spack_env/spack.yaml`.
You can add or remove packages from the `specs` list as needed or use different variants of these. 
For instance, if you want to use the debug variant (which builds Quandary and Petsc in debug mode) you can use `quandary@develop+test build_type=Debug`.
To build with `PetscInt`s set to 64-bit integers, use `quandary+int64`.
To use a specific version of Petsc instead of the latest release, you can do e.g. `quandary@develop^petsc@3.22.1`.
The `+test` variant (by default on in `.spack_env/spack.yaml`) adds python and pip to the Spack environment.
This allows `pip install` to install python packages to your Spack virtual environment rather than a global one.

## Building without Spack using CMake
If you don't want to use Spack to install all dependencies as explained above, you can follow these steps to install Petsc yourself, and then use CMake to install Quandary using the existing Petsc intallation, see below.

### Manually installing Petsc 
On MacOS, you can `brew install petsc`, or check out [https://petsc.org/release/] for the latest installation guide. As a quick start, you can try this: 
 * Download tarball for Petsc here [https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/] and extract `tar -xf petsc-<version>.tar.gz && cd petsc-<version>`
 * Configure Petsc with `./configure`. Check [https://petsc.org/release/install/install_tutorial] for optional arguments. For example, 
        `./configure --prefix=/YOUR/INSTALL/DIR --with-debugging=0 --with-fc=0 --with-cxx=mpicxx --with-cc=mpicc COPTFLAGS='-O3' CXXOPTFLAGS='-O3'`
 * The output of `./configure` reports on how to set the `PETSC_DIR` and `PETSC_ARCH` variables.
 You can export them or just note them, they are only needed to configure CMake within the PkgConfig step below.
      * `export PETSC_DIR=/YOUR/INSTALL/DIR`
      * `export PETSC_ARCH=/YOUR/ARCH/PREFIX`
 * Compile petsc with `make all check'


### Building with CMake
Given an existing Petsc installation, located in `$PETSC_DIR/$PETSC_ARCH`, you can build Quandary directly with CMake.
First tell PkgConfig (used by CMake) where to find your Petsc installation:
```
export PKG_CONFIG_PATH=$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig/:$PKG_CONFIG_PATH
```

Then build Quandary using:
```
mkdir build && cd build
cmake ..
make
```

To build in debug mode use `cmake -DCMAKE_BUILD_TYPE=Debug ..`.
Add the path to Quandary to your `PATH` variable with `export PATH=/path/to/quandary/:$PATH`, so your binary can be found.
Alternatively, you can install the Quandary executable in a specific path (such as the default `/usr/local/bin` to have it in your `PATH` automatically):
```
sudo cmake --install . --prefix /your/install/path
```

### Python dependencies and interface

Create a virtual environment (e.g. with conda, venv, ...) and then use `pip install -e .` to install the python dependencies and activate Quandary's python interface. For example, for Conda environments, do:
```
conda create --name myenv
conda activate myenv
pip install -e .
```
Or for venv, do:
```
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```
In additional to installing dependencies, the `pip install` command also enables that your python scripts can find Quandary's python interface functions, which are defined in `quandary.py`.


# Running
The C++ code builds into the executable `quandary`,
which takes one argument being the name of the test-case's configuration file. The file `config_template.cfg`, lists all possible configuration options and is filled with comments that should help users set up their own test case and match the options to the description in the user guide. 
* `quandary config_template.cfg` (serial execution)
* `mpirun -np 4 quandary config_template.cfg` (on 4 cores)

You can silence Quandary output by adding the `--quiet` argument to the above commands.

The `examples/` folder exemplifies the usage of Quandary's Python interface. 
* `python example_cnot.py`

# Building Documentation Locally

## Building Doxygen Locally
To build the Doxygen documentation locally:

1. Install Doxygen and Graphviz (e.g., `brew install doxygen graphviz` on macOS or `apt install doxygen graphviz` on Ubuntu)
2. Configure and build the docs:
   ```
   cd build
   cmake ..
   make quandary_doxygen
   ```
3. View the docs in your browser:
   ```
   open build/docs/doxygen/html/index.html
   ```

## Building User Guide Locally
To build the user guide documentation locally:

1. Install the required Python packages:
   ```
   pip install mkdocs pymdown-extensions mkdocs-bibtex mkdocs-macros-plugin mkdocs-jupyter
   ```
2. Build the MkDocs site:
   ```
   mkdocs build
   ```
3. To view the documentation in your browser while making changes (with live reload):
   ```
   mkdocs serve
   ```
   Then open http://127.0.0.1:8000/ in your browser.
4. The built site will be in the `site` directory.

# Tests


## Regression tests
Regression tests are defined in `tests/regression` and `tests/python` directories.

You can run all regression tests with:
```bash
pytest -m regression
```

Or run tests in a specific directory:
```bash
pytest tests/regression
pytest tests/python
```
See `tests/regression/README.md` for more information.

## Performance tests
Performance regression tests are defined in `tests/performance`.
The latest results from `main` are shown on this [performance dashboard](https://software.llnl.gov/quandary/dev/bench/).
See `tests/performance/README.md` for more information.

# Community and Contributing

Quandary is an open source project that is under heavy development. Contributions in all forms are very welcome, and can be anything from new features to bugfixes, documentation, or even discussions. Contributing is easy, work on your branch, create a pull request to `main` when you're good to go and the regression tests pass.

# Publications
* S. Guenther, N.A. Petersson, J.L. DuBois: "Quantum Optimal Control for Pure-State Preparation Using One Initial State", AVS Quantum Science, vol. 3, arXiv preprint <https://arxiv.org/abs/2106.09148> (2021)
* S. Guenther, N.A. Petersson, J.L. DuBois: "Quandary: An open-source C++ package for High-Performance Optimal Control of Open Quantum Systems", submitted to IEEE Supercomputing 2021, arXiv preprint <https://arxiv.org/abs/2110.10310> (2021)

# License

Quandary is distributed under the terms of the MIT license. All new contributions must be made under this license. See LICENSE, and NOTICE, for details. 

SPDX-License-Identifier: MIT
