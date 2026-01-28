"""Quandary Python interface for quantum optimal control."""

import warnings

# Import the new interface through the 'new' namespace
from . import new

# Import the legacy Quandary class and utility functions with deprecation warning
from .quandary import (
    Quandary,
    estimate_timesteps,
    eigen_and_reorder,
    get_resonances,
    lowering,
    number,
    map_to_oscillators,
    resolve_datadir,
    hamiltonians,
    plot_pulse,
    plot_expectedEnergy,
    plot_population,
    plot_results_1osc,
    timestep_richardson_est,
    execute,
    assemble_batch_script,
    infidelity_,
)

warnings.warn(
    "The legacy Quandary interface is deprecated and will be removed in a "
    "future version. Please migrate to the new interface: "
    "'from quandary.new import Quandary, QuandaryConfig'",
    DeprecationWarning,
    stacklevel=2
)

__all__ = [
    "Quandary",
    "estimate_timesteps",
    "eigen_and_reorder",
    "get_resonances",
    "lowering",
    "number",
    "map_to_oscillators",
    "resolve_datadir",
    "hamiltonians",
    "plot_pulse",
    "plot_expectedEnergy",
    "plot_population",
    "plot_results_1osc",
    "timestep_richardson_est",
    "execute",
    "assemble_batch_script",
    "infidelity_",
    "new",
]
