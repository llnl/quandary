# Quandary Parameter Reference

This document provides a comprehensive mapping between parameters in different Quandary interfaces.
Note that the CFG format is deprecated in favor of the TOML format and will be removed. See also [Config Options](config.md).

## System Parameters

| TOML Setting | Python Interface | CFG Format | Description |
|-------------|----------------|------------|-------------|
| `nlevels` | `Ne` + `Ng` | `nlevels` | Number of energy levels per oscillator |
| `nessential` | `Ne` | `nessential` | Number of essential levels per subsystem |
| `ntime` | Derived | `ntime` | Number of time steps used for time-integration |
| `dt` | Derived from `T`/`nsteps` | `dt` | Time step size (ns) |
| `transition_frequency` | `freq01` | `transfreq` | 01-transition frequencies [GHz] |
| `selfkerr` | `selfkerr` | `selfkerr` | Anharmonicities [GHz] |
| `crosskerr_coupling` | `crosskerr` | `crosskerr` | Cross-Kerr coupling strength [GHz] |
| `dipole_coupling` | `Jkl` | `Jkl` | Dipole-dipole coupling strength [GHz] |
| `rotation_frequency` | `rotfreq` | `rotfreq` | Rotational wave frequencies [GHz] |
| `decoherence` | `T1`/`T2` | `collapse_type` | Decoherence settings |
| `initial_condition` | `initialcondition` | `initialcondition` | Initial quantum state |
| `hamiltonian_file_Hsys` | `Hsys` | `hamiltonian_file_Hsys` | System Hamiltonian file |
| `hamiltonian_file_Hc` | `Hc_re`/`Hc_im` | `hamiltonian_file_Hc` | Control Hamiltonian file |

## Control Parameters

| TOML Setting | Python Interface | CFG Format | Description |
|-------------|----------------|------------|-------------|
| `parameterization` | `spline_order`, `nsplines` | `control_segments` | Control parameterization |
| `carrier_frequency` | `carrier_frequency` | `carrier_frequency` | Carrier frequencies |
| `initialization` | `initctrl_MHz`, `randomize_init_ctrl` | `control_init_type` | Control initialization |
| `amplitude_bound` | `maxctrl_MHz` | `control_bounds` | Control amplitude bounds |
| `zero_boundary_condition` | `control_enforce_BC` | `control_enforceBC` | Control boundary condition |

## Optimization Parameters

| TOML Setting | Python Interface | CFG Format | Description |
|-------------|----------------|------------|-------------|
| `target` | `targetgate`/`targetstate` | `optim_target` | Optimization target |
| `objective` | `costfunction` | `optim_objective` | Objective function |
| `weights` | N/A | `optim_weights` | Objective function weights |
| `tolerance` | `tol_infidelity`, `tol_costfunc`, etc. | `optim_tol_*` | Optimization tolerances |
| `maxiter` | `maxiter` | `optim_maxiter` | Maximum iterations |
| `tikhonov` | `gamma_tik0`, `gamma_tik0_interpolate` | `optim_regul` | Regularization parameters |
| `penalty` | `gamma_leakage`, `gamma_energy`, etc. | `optim_penalty_*` | Penalty parameters |

## Output Parameters

| TOML Setting | Python Interface | CFG Format | Description |
|-------------|----------------|------------|-------------|
| `directory` | `datadir` (method param) | `datadir` | Output directory |
| `observables` | N/A | `output` | Observables to output |
| `timestep_stride` | N/A | `output_frequency` | Output frequency |
| `optimization_stride` | `print_frequency_iter` | `optim_monitor_frequency` | Optimization output frequency |

## Solver Parameters

| TOML Setting | Python Interface | CFG Format | Description |
|-------------|----------------|------------|-------------|
| `runtype` | N/A | `runtype` | Run type |
| `usematfree` | `usematfree` | `usematfree` | Matrix-free solver setting |
| `linearsolver` | N/A | `linearsolver_type` | Linear solver configuration |
| `timestepper` | `timestepper` | `timestepper` | Time-stepping algorithm |
| `rand_seed` | `rand_seed` | `rand_seed` | Random seed |

## Defaults: Differences between Old and New Python Interface

The following defaults differ between the old `Quandary` dataclass and the new Python nanobind API.

| Old Python Name | Old Default | New Python Name | New Default | Notes |
|-----------------|-------------|-----------------|-------------|-------|
| `initctrl_MHz` | `10.0` MHz | `optimize(control_initialization_amplitude=)` | `0.01` GHz | Old: scaled by `1/1000/sqrt(2)/N_carriers`; New: passed directly |
| `control_enforce_BC` | `False` | `setup_quandary(control_zero_boundary_condition=)` | `True` | |
| `nsplines` | auto from `spline_knot_spacing=3.0` | `setup_quandary(control_parametrization_num=)` | `10` | |
| `tol_costfunc` | `1e-4` | `Setup.optim_tol_final_cost` | `1e-8` | |
| `gamma_leakage` | `0.1` | `Setup.optim_penalty_leakage` | `0.0` | |
| `gamma_energy` | `0.1` | `Setup.optim_penalty_energy` | `0.0` | |
| `gamma_dpdm` | `0.01` | `Setup.optim_penalty_dpdm` | `0.0` | |
| `print_frequency_iter` | `1` | `Setup.output_optimization_stride` | `10` | |
| n/a (hardcoded in `__dump`) | `20` | `Setup.linearsolver_maxiter` | `10` | |
| `rand_seed` | `None` (random) | `Setup.rand_seed` | `1` (deterministic) | Only differs when not explicitly set |
