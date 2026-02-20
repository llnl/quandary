"""New Quandary Python interface for quantum optimal control using nanobind."""

# Initialize MPI via mpi4py before loading the C++ extension.
# This ensures MPI is initialized once and managed by Python's lifecycle.
# MPI_Init called twice is undefined behavior, so this must happen before any C++ code.
import mpi4py.MPI  # noqa: F401, E402

# Re-export the nanobind implementation from the parent package
from .._quandary_impl import (
    # Exceptions
    ValidationError as ValidationError,
    # Configuration
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
)

# Python subclass of Setup with __repr__ and improved TypeError messages
from ._structs import Setup

from .._quandary_impl import (
    InitialConditionSettings,
    OptimTargetSettings,
    ControlParameterizationSettings,
    ControlInitializationSettings,
)


# Functional API (primary interface)
from .runner import optimize, simulate, evaluate_controls  # noqa: F401, E402

# Results
from .results import Results, get_results  # noqa: F401, E402

# Quantum operators and Hamiltonians
from .quantum_operators import (  # noqa: F401, E402
    lowering,
    number,
    map_to_oscillators,
    hamiltonians,
    get_resonances,
)

# Time estimation utilities
from .time_estimation import (  # noqa: F401, E402
    estimate_timesteps,
    timestep_richardson_est,
)

# Visualization utilities
from .visualization import (  # noqa: F401, E402
    plot_pulse,
    plot_expectedEnergy,
    plot_population,
    plot_results_1osc,
)

# General utilities
from .utils import (  # noqa: F401, E402
    gate_infidelity,
    state_infidelity,
    downsample_pulses,
)

# Setup helpers (factory/configuration functions)
from .setup_helpers import (  # noqa: F401, E402
    setup_quandary,
    resolve_output_dir,
)

# Define public API
__all__ = [
    # Exceptions
    "ValidationError",
    # Configuration classes and types
    "Setup",
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
    "optimize",
    "simulate",
    "evaluate_controls",
    # Results
    "Results",
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
    "gate_infidelity",
    "state_infidelity",
    "downsample_pulses",
    # Setup helpers
    "setup_quandary",
    "resolve_output_dir",
]
