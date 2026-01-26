"""Quandary Python interface for quantum optimal control."""

from __future__ import annotations

from .quandary_ext import (
  # Main API
  ConfigInput as _ConfigInput,
  Config,
  run,
  run_from_file,
  # Enums
  DecoherenceType,
  InitialConditionType,
  TargetType,
  ObjectiveType,
  LinearSolverType,
  RunType,
  ControlType,
  ControlInitializationType,
  TimeStepperType,
  GateType,
  OutputType,
  # Structs
  InitialConditionSettings,
  OptimTargetSettings,
  ControlParameterizationSettings,
  ControlInitializationSettings,
)


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

__all__ = [
  # Main API
  "Quandary",
  "Config",
  "run",
  "run_from_file",
  # Enums
  "DecoherenceType",
  "InitialConditionType",
  "TargetType",
  "ObjectiveType",
  "LinearSolverType",
  "RunType",
  "ControlType",
  "ControlInitializationType",
  "TimeStepperType",
  "GateType",
  "OutputType",
  # Structs
  "InitialConditionSettings",
  "OptimTargetSettings",
  "ControlParameterizationSettings",
  "ControlInitializationSettings",
]
