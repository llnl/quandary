# Quandary Parameter Reference

This document provides a comprehensive mapping between parameters in different Quandary interfaces.
Note that the CFG format is deprecated in favor of the TOML format and will be removed. See also [Config Options](../config.md).

## System Parameters

| TOML Setting | Python Interface | CFG Format | Description |
|-------------|----------------|------------|-------------|
| `nlevels` | `Ne` + `Ng` | `nlevels` | Number of energy levels per oscillator |
| `nessential` | `Ne` | N/A | Number of essential levels per subsystem |
| `ntime` | Derived | `ntime` | Number of time steps used for time-integration |
| `dt` | Derived from `T`/`nsteps` | `dt` | Time step size (ns) |
| `transition_frequency` | `freq01` | `transfreq` | 01-transition frequencies [GHz] |
| `selfkerr` | `selfkerr` | `selfkerr` | Anharmonicities [GHz] |
| `crosskerr` | `crosskerr` | `crosskerr` | ZZ coupling strength [GHz] |
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
| `zero_boundary_condition` | `control_enforce_BC` | `control_zero_boundary_condition` | Control boundary condition |

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
| `directory` | N/A | `datadir` | Output directory |
| `observables` | N/A | `output` | Observables to output |
| `timestep_stride` | N/A | `output_timestep_stride` | Output frequency |
| `optimization_stride` | `print_frequency_iter` | `output_optimization_stride` | Optimization output frequency |

## Solver Parameters

| TOML Setting | Python Interface | CFG Format | Description |
|-------------|----------------|------------|-------------|
| `runtype` | N/A | `runtype` | Run type |
| `usematfree` | `usematfree` | `usematfree` | Matrix-free solver setting |
| `linearsolver` | N/A | `linearsolver_type` | Linear solver configuration |
| `timestepper` | `timestepper` | `timestepper` | Time-stepping algorithm |
| `rand_seed` | `rand_seed` | `rand_seed` | Random seed |
