"""Quandary runner class."""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from typing import TYPE_CHECKING, Optional

from .. import _quandary_impl
from .._quandary_impl import QuandaryConfig, Config, DecoherenceType
from .results import get_results as _get_results

if TYPE_CHECKING:
    from .results import QuandaryResults

logger = logging.getLogger(__name__)


class Quandary:
    """Runner for Quandary simulations and optimizations.

    Takes a QuandaryConfig, validates it, and provides methods to run
    simulations and retrieve results.

    The configuration is validated at construction time. After construction,
    the Quandary object is immutable.

    Example:
        config = QuandaryConfig()
        config.nlevels = [2]
        config.ntime = 1000
        config.dt = 0.01
        config.transfreq = [4.1]
        config.runtype = RunType.SIMULATION

        quandary = Quandary(config)
        quandary.run()
        results = quandary.get_results()
    """

    def __init__(self, config: QuandaryConfig, quiet: bool = False):
        """Create a Quandary runner from a configuration.

        Args:
            config: A QuandaryConfig object with all required fields set.
            quiet: If True, suppress console output during validation.

        Raises:
            RuntimeError: If the configuration is invalid.
        """
        self._input = config
        self._config = Config(config, quiet)

    @property
    def config(self) -> Config:
        """Get the validated Config object."""
        return self._config

    @property
    def n_initial_conditions(self) -> int:
        """Number of initial conditions."""
        return self._config.n_initial_conditions

    @property
    def output_directory(self) -> str:
        """Output directory path."""
        return self._config.output_directory

    @property
    def decoherence_type(self):
        """Decoherence type (NONE, DECAY, DEPHASE, or BOTH)."""
        return self._config.decoherence_type

    def to_toml(self) -> str:
        """Serialize the configuration to TOML format.

        Returns:
            String containing the configuration in TOML format.
        """
        return self._config.to_toml()

    def run(self, quiet: bool = False) -> int:
        """Run the simulation or optimization.

        Args:
            quiet: If True, suppress console output.

        Returns:
            0 on success, non-zero on error.
        """
        return _quandary_impl.run(self._config, quiet)

    def run_mpi(
        self,
        n_procs: int,
        quiet: bool = False,
        mpi_exec: str = "mpirun",
        python_exec: Optional[str] = None,
        working_dir: Optional[str] = None,
    ) -> "QuandaryResults":
        """Run the simulation or optimization using MPI.

        This method writes the configuration to a TOML file, then runs Quandary
        via Python as a subprocess with MPI. This is useful for Jupyter notebooks
        where the Python process cannot directly use MPI.

        Args:
            n_procs: Number of MPI processes to use.
            quiet: If True, suppress console output.
            mpi_exec: MPI launcher command (e.g., "mpirun", "srun").
            python_exec: Path to the Python executable. Defaults to sys.executable.
            working_dir: Working directory for the subprocess. Defaults to
                the output_directory.

        Returns:
            QuandaryResults containing all parsed output data.

        Raises:
            subprocess.CalledProcessError: If the Quandary process fails.

        Example:
            >>> quandary = Quandary(config)
            >>> results = quandary.run_mpi(n_procs=4)
            >>> print(f"Infidelity: {results.infidelity}")
        """
        # Use output_directory as working directory if not specified
        if working_dir is None:
            working_dir = self.output_directory

        # Use current Python interpreter if not specified
        if python_exec is None:
            python_exec = sys.executable

        # Create the working directory if it doesn't exist
        os.makedirs(working_dir, exist_ok=True)

        # Write the TOML config file (use absolute path for subprocess)
        config_file = os.path.abspath(os.path.join(working_dir, "config.toml"))
        toml_content = self.to_toml()
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

        # Return the results
        return self.get_results()

    def get_results(self) -> "QuandaryResults":
        """Load results from the output directory.

        Returns:
            QuandaryResults containing all parsed output data.

        Example:
            >>> quandary = Quandary(config)
            >>> quandary.run()
            >>> results = quandary.get_results()
            >>> print(f"Infidelity: {results.infidelity}")
        """
        return _get_results(
            datadir=self.output_directory,
            lindblad=self.decoherence_type != DecoherenceType.NONE,
            n_init=self.n_initial_conditions,
        )
