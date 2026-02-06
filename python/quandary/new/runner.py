"""Functional API for running Quandary simulations and optimizations.

MPI Lifecycle Management:
    MPI initialization and finalization are managed by mpi4py, which is imported
    in __init__.py. This allows multiple run() calls in the same Python process
    (e.g., Jupyter notebooks, parameter sweeps) without MPI re-initialization errors.

    The C++ runQuandary function is called with finalize=false to prevent it from
    calling MPI_Finalize, allowing mpi4py to manage finalization at process exit.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from typing import TYPE_CHECKING, Optional

from .. import _quandary_impl
from .._quandary_impl import Config, QuandaryConfig
from .results import get_results as _get_results, QuandaryResults

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


def run(config: QuandaryConfig, quiet: bool = False) -> QuandaryResults:
    """Run a Quandary simulation or optimization.

    This function validates the configuration, runs the simulation/optimization,
    and returns results with the validated configuration attached.

    If already running in an MPI context (launched with mpirun), MPI will be
    used automatically. Otherwise runs serially.

    Args:
        config: A QuandaryConfig object with all required fields set.
        quiet: If True, suppress console output.

    Returns:
        QuandaryResults containing output data and the validated configuration.

    Raises:
        RuntimeError: If the configuration is invalid or execution fails.

    Example:
        >>> config = QuandaryConfig()
        >>> config.nlevels = [2]
        >>> config.ntime = 1000
        >>> config.dt = 0.01
        >>> config.transition_frequency = [4.1]
        >>> config.runtype = RunType.SIMULATION
        >>> results = run(config)
        >>> print(f"Infidelity: {results.infidelity}")
        >>> print(results.config.to_toml())  # See validated config
    """
    # Validate configuration
    validated_config = Config(config, quiet)

    # Run simulation/optimization
    return_code = _quandary_impl.run(validated_config, quiet)

    if return_code != 0:
        raise RuntimeError(f"Quandary execution failed with return code {return_code}")

    # Load results with validated config
    results = _get_results(validated_config)

    return results


def run_mpi(
    config: QuandaryConfig,
    n_procs: int,
    quiet: bool = False,
    mpi_exec: str = "mpirun",
    python_exec: Optional[str] = None,
    working_dir: str = ".",
) -> QuandaryResults:
    """Run Quandary with MPI via subprocess.

    This function is designed for Jupyter notebooks and other contexts where
    the Python process cannot directly use MPI. It writes the configuration
    to a TOML file, spawns a subprocess with the MPI launcher, and returns
    the results.

    For regular Python scripts, use mpirun to launch your script and call
    run() directly instead.

    Args:
        config: A QuandaryConfig object with all required fields set.
        n_procs: Number of MPI processes to use.
        quiet: If True, suppress console output.
        mpi_exec: MPI launcher command (e.g., "mpirun", "srun"). Default: "mpirun".
        python_exec: Path to Python executable. Defaults to sys.executable.
        working_dir: Working directory for subprocess. Defaults to current directory.

    Returns:
        QuandaryResults containing output data and the validated configuration.

    Raises:
        subprocess.CalledProcessError: If the Quandary process fails.

    Example:
        >>> # In Jupyter notebook
        >>> config = create_optimization_config(nessential=[2], transition_frequency=[4.1], final_time=50.0)
        >>> results = run_mpi(config, n_procs=4)
        >>> print(f"Infidelity: {results.infidelity}")
    """
    # Validate configuration first
    validated_config = Config(config, quiet)

    # Use current Python interpreter if not specified
    if python_exec is None:
        python_exec = sys.executable

    # Ensure working_dir is a string
    if not isinstance(working_dir, str):
        raise ValueError(f"working_dir must be a string, got {type(working_dir)}")

    # Create the output directory if it doesn't exist
    os.makedirs(validated_config.output_directory, exist_ok=True)

    # Write the TOML config file to the output directory (use absolute path for subprocess)
    config_file = os.path.abspath(os.path.join(validated_config.output_directory, "config.toml"))
    toml_content = validated_config.to_toml()
    with open(config_file, "w") as f:
        f.write(toml_content)

    # Python code to run Quandary from the TOML file
    python_code = f'from quandary.new import run_from_file; run_from_file("{config_file}", quiet={quiet})'

    # Build the command
    cmd = [mpi_exec, "-n", str(n_procs), python_exec, "-c", python_code]

    logger.info(f"Running MPI with {n_procs} processes")

    # Run the subprocess
    result = subprocess.run(
        cmd,
        cwd=working_dir,
        capture_output=quiet,
        text=True,
    )

    if result.returncode != 0:
        if quiet and result.stderr:
            logger.error(f"Quandary failed: {result.stderr}")
        result.check_returncode()

    # Load results with validated config
    results = _get_results(validated_config)

    return results
