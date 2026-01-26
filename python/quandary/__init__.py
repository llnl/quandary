"""Quandary Python interface for quantum optimal control."""

from .quandary_ext import (
  # Main API
  ConfigInput as _ConfigInput,
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
  """Configuration input for a Quandary simulation or optimization.

  All fields are optional. Set only the fields you need.

  Example:
    config = ConfigInput(
      nlevels=[2],
      ntime=1000,
      dt=0.01,
      transfreq=[4.1],
      runtype=RunType.SIMULATION,
    )
    run(config)
  """

  def __init__(self, **kwargs):
    super().__init__()
    for key, value in kwargs.items():
      if not hasattr(self, key):
        raise AttributeError(f"ConfigInput has no field '{key}'")
      setattr(self, key, value)

__all__ = [
  # Main API
  "Quandary",
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
