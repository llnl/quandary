"""New Quandary Python interface for quantum optimal control using nanobind."""

# Initialize MPI and PETSc via mpi4py/petsc4py before any C++ calls
# This ensures both are initialized once and managed by Python's lifecycle
import mpi4py.MPI  # noqa: F401
import petsc4py

# Re-export the nanobind implementation from the parent package
from .._quandary_impl import (
    # Configuration
    QuandaryConfig,
    Config,
    # Run function from C++ (used internally)
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

# Functional API (primary interface)
from .runner import run, validate  # noqa: F401

# Results
from .results import QuandaryResults, get_results  # noqa: F401

# Quantum operators and Hamiltonians
from .quantum_operators import (  # noqa: F401
    lowering,
    number,
    map_to_oscillators,
    hamiltonians,
    get_resonances,
)

# Time estimation utilities
from .time_estimation import (  # noqa: F401
    estimate_timesteps,
    timestep_richardson_est,
)

# Visualization utilities
from .visualization import (  # noqa: F401
    plot_pulse,
    plot_expectedEnergy,
    plot_population,
    plot_results_1osc,
)

# General utilities
from .utils import (  # noqa: F401
    eval_controls,
    infidelity_,
    downsample_pulses,
)

# Config builders (factory functions)
from .config_builders import (
    create_simulation_config,
    create_optimization_config,
)

# Define public API
__all__ = [
    # Configuration classes and types
    "QuandaryConfig",
    "Config",
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
    # Settings structs
    "InitialConditionSettings",
    "OptimTargetSettings",
    "ControlParameterizationSettings",
    "ControlInitializationSettings",
    # Runner functions
    "run",
    "validate",
    # Results
    "QuandaryResults",
    "get_results",
    # Quantum operators
    "lowering",
    "number",
    "map_to_oscillators",
    "hamiltonians",
    "get_resonances",
    # Time estimation
    "estimate_timesteps",
    "timestep_richardson_est",
    # Visualization
    "plot_pulse",
    "plot_expectedEnergy",
    "plot_population",
    "plot_results_1osc",
    # Utilities
    "eval_controls",
    "infidelity_",
    "downsample_pulses",
    # Config builders
    "create_simulation_config",
    "create_optimization_config",
]

petsc4py.init()
