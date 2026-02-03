"""New Quandary Python interface for quantum optimal control using nanobind."""

# Re-export the nanobind implementation from the parent package
from .._quandary_impl import (  # noqa: F401
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
from .quandary import Quandary  # noqa: F401
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
    infidelity_,
    downsample_pulses,
)
