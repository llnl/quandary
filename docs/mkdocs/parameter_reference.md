# Quandary Parameter Reference

This document provides a comprehensive mapping between parameters in the TOML configuration and the Python interface. See also [Config Options](config.md).

## System Parameters

| TOML Setting | Python Interface | Description |
|-------------|----------------|-------------|
| `nlevels` | `Ne` + `Ng` | Number of energy levels per oscillator |
| `nessential` | `Ne` | Number of essential levels per subsystem |
| `ntime` | Derived | Number of time steps used for time-integration |
| `dt` | Derived from `T`/`nsteps` | Time step size (ns) |
| `transition_frequency` | `freq01` | 01-transition frequencies [GHz] |
| `selfkerr` | `selfkerr` | Anharmonicities [GHz] |
| `crosskerr_coupling` | `crosskerr` | Cross-Kerr coupling strength [GHz] |
| `dipole_coupling` | `Jkl` | Dipole-dipole coupling strength [GHz] |
| `rotation_frequency` | `rotfreq` | Rotational wave frequencies [GHz] |
| `decoherence` | `T1`/`T2` | Decoherence settings |
| `initial_condition` | `initialcondition` | Initial quantum state |
| `hamiltonian_file_Hsys` | `Hsys` | System Hamiltonian file |
| `hamiltonian_file_Hc` | `Hc_re`/`Hc_im` | Control Hamiltonian file |

## Control Parameters

| TOML Setting | Python Interface | Description |
|-------------|----------------|-------------|
| `parameterization` | `spline_order`, `nsplines` | Control parameterization |
| `carrier_frequency` | `carrier_frequency` | Carrier frequencies |
| `initialization` | `initctrl_MHz`, `randomize_init_ctrl` | Control initialization |
| `amplitude_bound` | `maxctrl_MHz` | Control amplitude bounds |
| `zero_boundary_condition` | `control_enforce_BC` | Control boundary condition |

## Optimization Parameters

| TOML Setting | Python Interface | Description |
|-------------|----------------|-------------|
| `target` | `targetgate`/`targetstate` | Optimization target |
| `objective` | `costfunction` | Objective function |
| `weights` | N/A | Objective function weights |
| `tolerance` | `tol_infidelity`, `tol_costfunc`, etc. | Optimization tolerances |
| `maxiter` | `maxiter` | Maximum iterations |
| `tikhonov` | `gamma_tik0`, `gamma_tik0_interpolate` | Regularization parameters |
| `penalty` | `gamma_leakage`, `gamma_energy`, etc. | Penalty parameters |

## Output Parameters

| TOML Setting | Python Interface | Description |
|-------------|----------------|-------------|
| `directory` | `datadir` (method param) | Output directory |
| `observables` | N/A | Observables to output |
| `timestep_stride` | N/A | Output frequency |
| `optimization_stride` | `print_frequency_iter` | Optimization output frequency |

## Solver Parameters

| TOML Setting | Python Interface | Description |
|-------------|----------------|-------------|
| `runtype` | N/A | Run type |
| `usematfree` | `usematfree` | Matrix-free solver setting |
| `linearsolver` | N/A | Linear solver configuration |
| `timestepper` | `timestepper` | Time-stepping algorithm |
| `rand_seed` | `rand_seed` | Random seed |
