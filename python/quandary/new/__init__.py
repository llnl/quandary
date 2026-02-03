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

# Quantum operators and Hamiltonians
from .quantum_operators import (
    lowering,
    number,
    map_to_oscillators,
    hamiltonians,
    get_resonances,
)

# Time estimation utilities
from .time_estimation import (
    estimate_timesteps,
    timestep_richardson_est,
)

# Visualization utilities
from .visualization import (
    plot_pulse,
    plot_expectedEnergy,
    plot_population,
    plot_results_1osc,
)

# General utilities
from .utils import (
    infidelity_,
    downsample_pulses,
)
