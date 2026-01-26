"""Quandary Python interface for quantum optimal control."""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

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


@dataclass
class QuandaryResults:
  """Results from a Quandary simulation or optimization.

  Attributes:
    time: Array of time points [ns].
    pt: Control pulses p(t) for each oscillator [MHz]. Access: pt[oscillator][time_index].
    qt: Control pulses q(t) for each oscillator [MHz]. Access: qt[oscillator][time_index].
    ft: Lab-frame control pulses f(t) for each oscillator [MHz]. Access: ft[oscillator][time_index].
    params: Control parameters (B-spline coefficients).
    infidelity: Final infidelity (1 - fidelity).
    optim_history: Optimization history with keys: 'iter', 'objective', 'gradient',
        'ls_step', 'fidelity', 'cost', 'tikhonov', 'penalty', 'state_variation',
        'energy', 'control_variation'.
    expected_energy: Expected energy evolution for each oscillator and initial condition.
        Access: expected_energy[oscillator][initial_condition][time_index].
    population: Population evolution for each oscillator and initial condition.
        Access: population[oscillator][initial_condition][level, time_index].
    rho_final: Final density matrix/state for each initial condition.
        Access: rho_final[:, initial_condition].
  """
  time: np.ndarray = field(default_factory=lambda: np.array([]))
  pt: List[np.ndarray] = field(default_factory=list)
  qt: List[np.ndarray] = field(default_factory=list)
  ft: List[np.ndarray] = field(default_factory=list)
  params: np.ndarray = field(default_factory=lambda: np.array([]))
  infidelity: float = 1.0
  optim_history: Dict[str, np.ndarray] = field(default_factory=dict)
  expected_energy: List[List[np.ndarray]] = field(default_factory=list)
  population: List[List[np.ndarray]] = field(default_factory=list)
  rho_final: np.ndarray = field(default_factory=lambda: np.array([]))


def get_results(
    datadir: str = "./",
    config: Optional[Config] = None,
) -> QuandaryResults:
  """Load results from Quandary output files.

  This function parses output files from a Quandary run and returns them
  in a structured format. If a Config is provided, its computed values
  (n_initial_conditions, lindblad) are used; otherwise they are auto-detected.

  Args:
    datadir: Directory containing Quandary output files.
    config: Optional Config object to get n_initial_conditions and lindblad from.

  Returns:
    QuandaryResults containing all parsed output data.

  Example:
    >>> config = Config(Quandary(nlevels=[2], ntime=100, dt=0.01))
    >>> run(config)
    >>> results = get_results(config.output_directory, config)
    >>> print(f"Infidelity: {results.infidelity}")
  """
  results = QuandaryResults()

  # Auto-detect from files (always needed for rho_files)
  control_files = sorted(glob.glob(os.path.join(datadir, "control*.dat")))
  rho_files = sorted(glob.glob(os.path.join(datadir, "rho_Re.iinit*.dat")))
  n_osc = len(control_files)

  # Get lindblad and n_init from Config if provided, otherwise auto-detect
  if config is not None:
    lindblad = config.lindblad
    n_init = config.n_initial_conditions
  else:
    n_init = len(rho_files)
    # Auto-detect whether Lindblad solver was used from config_log.toml
    lindblad = False
    config_log = os.path.join(datadir, "config_log.toml")
    if os.path.exists(config_log):
      try:
        with open(config_log) as f:
          for line in f:
            if "type" in line and any(x in line for x in ['"decay"', '"dephase"', '"both"']):
              lindblad = True
              break
      except Exception:
        pass

  # For Lindblad, we only want diagonal initial conditions for some outputs
  n_init_diag = n_init if not lindblad else int(np.sqrt(n_init))

  # Read control parameters (params.dat)
  params_file = os.path.join(datadir, "params.dat")
  if os.path.exists(params_file):
    try:
      results.params = np.loadtxt(params_file)
    except Exception:
      pass

  # Read optimization history (optim_history.dat)
  optim_file = os.path.join(datadir, "optim_history.dat")
  if os.path.exists(optim_file):
    try:
      data = np.loadtxt(optim_file)
      if data.ndim == 1:
        data = data.reshape(1, -1)
      results.optim_history = {
        "iter": data[:, 0],
        "objective": data[:, 1],
        "gradient": data[:, 2],
        "ls_step": data[:, 3],
        "fidelity": data[:, 4],
        "cost": data[:, 5],
        "tikhonov": data[:, 6],
        "penalty": data[:, 7],
        "state_variation": data[:, 8],
        "energy": data[:, 9],
        "control_variation": data[:, 10] if data.shape[1] > 10 else np.zeros_like(data[:, 0]),
      }
      # Infidelity from last iteration
      results.infidelity = 1.0 - data[-1, 4]
    except Exception:
      pass

  # Read control pulses for each oscillator
  # Convert from GHz (file) to MHz (output)
  ghz_to_mhz = 1e3
  for iosc in range(n_osc):
    ctrl_file = os.path.join(datadir, f"control{iosc}.dat")
    if os.path.exists(ctrl_file):
      try:
        data = np.loadtxt(ctrl_file)
        if iosc == 0:
          results.time = data[:, 0]
        results.pt.append(data[:, 1] * ghz_to_mhz)
        results.qt.append(data[:, 2] * ghz_to_mhz)
        results.ft.append(data[:, 3] * ghz_to_mhz)
      except Exception:
        pass

  # Read expected energy for each oscillator and initial condition
  for iosc in range(n_osc):
    osc_expected = []
    for iinit in range(n_init_diag):
      # For Lindblad, only read diagonal initial conditions
      iid = iinit if not lindblad else iinit * n_init_diag + iinit
      filename = os.path.join(datadir, f"expected{iosc}.iinit{iid:04d}.dat")
      if os.path.exists(filename):
        try:
          data = np.loadtxt(filename)
          osc_expected.append(data[:, 1])
        except Exception:
          pass
    if osc_expected:
      results.expected_energy.append(osc_expected)

  # Read population for each oscillator and initial condition
  for iosc in range(n_osc):
    osc_population = []
    for iinit in range(n_init_diag):
      iid = iinit if not lindblad else iinit * n_init_diag + iinit
      filename = os.path.join(datadir, f"population{iosc}.iinit{iid:04d}.dat")
      if os.path.exists(filename):
        try:
          data = np.loadtxt(filename)
          # Population data: first column is time, rest are level populations
          osc_population.append(data[:, 1:].T)
        except Exception:
          pass
    if osc_population:
      results.population.append(osc_population)

  # Read final state/density matrix for each initial condition
  # First, determine dimensions from the first rho file
  if rho_files and n_init > 0:
    try:
      # Read first file to get dimensions
      first_rho = np.loadtxt(rho_files[0], skiprows=1)
      ndim = first_rho.shape[1] - 1  # First column is time

      # Initialize rho_final array
      results.rho_final = np.zeros((ndim, n_init), dtype=complex)

      # Read each initial condition
      for iinit in range(n_init):
        file_index = f"{iinit:04d}"
        re_file = os.path.join(datadir, f"rho_Re.iinit{file_index}.dat")
        im_file = os.path.join(datadir, f"rho_Im.iinit{file_index}.dat")

        if os.path.exists(re_file):
          try:
            re_data = np.loadtxt(re_file, skiprows=1)
            # Take last time step, skip time column
            results.rho_final[:, iinit] = re_data[-1, 1:]
          except Exception:
            pass

        if os.path.exists(im_file):
          try:
            im_data = np.loadtxt(im_file, skiprows=1)
            results.rho_final[:, iinit] += 1j * im_data[-1, 1:]
          except Exception:
            pass
    except Exception:
      pass

  return results

__all__ = [
  # Main API
  "Quandary",
  "Config",
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
