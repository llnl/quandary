"""New Quandary Python interface for quantum optimal control using nanobind."""

# Re-export the nanobind implementation from the parent package
from .._quandary_impl import (
    # Configuration
    QuandaryConfig,
    # Run functions
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

# Import our implementation
from .quandary import Quandary
from .results import QuandaryResults, get_results

__all__ = [
  # Main API
  "QuandaryConfig",
  "Quandary",
  "QuandaryResults",
  "run",
  "run_from_file",
  "get_results",
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
