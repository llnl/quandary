"""Quandary configuration and main interface class."""

from __future__ import annotations

import logging
import os
import subprocess
from typing import TYPE_CHECKING, Optional

from .quandary_ext import ConfigInput as _ConfigInput, Config

if TYPE_CHECKING:
  from .results import QuandaryResults

logger = logging.getLogger(__name__)


class Quandary(_ConfigInput):
  """Configuration for a Quandary simulation or optimization.

  Accepts all configuration fields as keyword arguments. The configuration
  is validated immediately at construction. Computed properties
  (n_initial_conditions, lindblad, output_directory, etc.) are available
  after construction via automatic forwarding to the validated Config.

  The object is immutable after construction - attempting to change fields
  will raise an error.

  Example:
    config = Quandary(
      nlevels=[2],
      ntime=1000,
      dt=0.01,
      transfreq=[4.1],
      runtype=RunType.SIMULATION,
    )
    print(config.n_initial_conditions)  # Computed value available
    run(config)
  """

  def __init__(self, quiet: bool = False, **kwargs):
    # Use object.__setattr__ to bypass our __setattr__ during init
    object.__setattr__(self, '_frozen', False)

    super().__init__()

    for key, value in kwargs.items():
      if not hasattr(self, key):
        raise AttributeError(f"Quandary has no field '{key}'")
      setattr(self, key, value)

    # Eagerly create and validate Config
    object.__setattr__(self, '_config', Config(self, quiet))

    # Freeze the object - no more changes allowed
    object.__setattr__(self, '_frozen', True)

  @property
  def config(self) -> Config:
    """Get validated Config with computed values."""
    return self._config

  def __setattr__(self, name: str, value):
    """Prevent modification after construction."""
    if getattr(self, '_frozen', False) and not name.startswith('_'):
      raise AttributeError(
        f"Quandary is immutable after construction. "
        f"Cannot modify '{name}'. Create a new Quandary instead."
      )
    super().__setattr__(name, value)

  def __getattr__(self, name: str):
    """Forward unknown attributes to the validated Config."""
    # Avoid infinite recursion for private attributes
    if name.startswith('_'):
      raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")
    # Forward to Config
    return getattr(self._config, name)

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
    from .quandary_ext import run as _run
    return _run(self._config, quiet)

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
      >>> config = Quandary(nlevels=[2], ntime=100, dt=0.01, ...)
      >>> results = config.run_mpi(n_procs=4)
      >>> print(f"Infidelity: {results.infidelity}")
    """
    import sys

    # Use output_directory as working directory if not specified
    if working_dir is None:
      working_dir = str(self.output_directory)

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
    python_code = f'from quandary import run_from_file; run_from_file("{config_file}", quiet={quiet})'

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
      >>> config = Quandary(nlevels=[2], ntime=100, dt=0.01, ...)
      >>> run(config)
      >>> results = config.get_results()
      >>> print(f"Infidelity: {results.infidelity}")
    """
    from .results import get_results as _get_results
    return _get_results(
      datadir=self.output_directory,
      lindblad=self.lindblad,
      n_init=self.n_initial_conditions,
    )
