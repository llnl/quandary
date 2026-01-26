"""Quandary configuration and main interface class."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .quandary_ext import ConfigInput as _ConfigInput, Config

if TYPE_CHECKING:
  from .results import QuandaryResults


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
