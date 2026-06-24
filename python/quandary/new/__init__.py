"""New Quandary Python interface for quantum optimal control using nanobind."""

# Initialize MPI via mpi4py before loading the C++ extension.
# This ensures MPI is initialized once and managed by Python's lifecycle.
import mpi4py.MPI  # noqa: F401, E402
# Register PETSc finalization to run before mpi4py's MPI_Finalize at exit.
# atexit handlers fire in reverse registration order, so registering after
# mpi4py ensures PETSc is finalized while MPI is still active.
import atexit
from .._quandary_impl import _finalize_petsc
atexit.register(_finalize_petsc)

# Re-export the nanobind implementation from the parent package
from .._quandary_impl import (
    ValidationError as ValidationError,
    run_from_file,
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

# Types
from .types import (  # noqa: F401, E402
    Config,
    InitialConditionSettings,
    OptimTargetSettings,
    ControlParameterizationSettings,
    ControlInitializationSettings,
)

# Configuration functions
from .config import (  # noqa: F401, E402
    create_config,
    load_config,
    resolve_output_dir,
    set_target,
    set_initial_condition,
    set_controls,
    set_decoherence,
)

# Runner functions
from .run import (  # noqa: F401, E402
    optimize,
    simulate,
    evaluate_controls,
)

# Results and visualization
from .results import (  # noqa: F401, E402
    Results,
    get_results,
    plot_pulse,
    plot_expectedEnergy,
    plot_population,
    plot_results_1osc,
)

# Quantum operators, Hamiltonians, and utilities
from .physics import (  # noqa: F401, E402
    lowering,
    number,
    hamiltonians,
    get_resonances,
    gate_infidelity,
    state_infidelity,
    estimate_timestep_size,
    timestep_richardson_est,
    fit_bspline0,
    fit_bspline2nd,
)

# Define public API
__all__ = [
    # Exceptions
    "ValidationError",
    # Configuration classes and types
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
    # Configuration 
    "create_config",
    "load_config",
    "resolve_output_dir",
    "set_target",
    "set_controls",
    "set_initial_condition",
    "set_decoherence",
    # Runner functions
    "optimize",
    "simulate",
    "evaluate_controls",
    # Results
    "Results",
    "get_results",
    "plot_pulse",
    "plot_expectedEnergy",
    "plot_population",
    "plot_results_1osc",
    # Quantum operators and utilities
    "lowering",
    "number",
    "hamiltonians",
    "get_resonances",
    "gate_infidelity",
    "state_infidelity",
    "fit_bspline0",
    "fit_bspline2nd",
    "estimate_timestep_size",
    "timestep_richardson_est",
]
