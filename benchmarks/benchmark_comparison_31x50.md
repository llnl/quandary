# PETSc Benchmark Comparison: System 31×50 (Lindblad, N=2,402,500)

## Run Summary
Note: 4 GPUs means 4 APUs and 4 ranks (-n 4).

| Timer | 4 CPU ranks | 62 CPU ranks | 4 GPU ranks (4 MI300A) |
|---|---:|---:|---:|
| **Quandary Used Time** | **563.9s** | **119.0s** | **58.6s** |
| PETSc Wall clock | 591.9s | 175.3s | 74.5s |
| PETSc KSPSolve | 521.0s | 114.2s | 53.4s |
| PETSc MatMult | 435.8s | 94.7s | ~23.4s (est.) |
| PETSc Mflop/s (KSPSolve) | 7,212 | 42,385 | 74,314 |

## Speedups

| Timer | 4 CPU → 62 CPU | 4 CPU → 4 GPU | 62 CPU → 4 GPU |
|---|---:|---:|---:|
| **Quandary Used Time** | 4.7× | 9.6× | 2.0× |
| PETSc Wall clock | 3.4× | 7.9× | 2.4× |
| PETSc KSPSolve | 4.6× | 9.8× | 2.1× |
| PETSc MatMult | 4.6× | ~18.6× | ~4.0× |

## Detailed Event Breakdown

Times in seconds. GPU times estimated from %T × Used Time (58.6s) since PETSc reports `n/a` without `-log_view_gpu_time`.

| Event | 4 CPU | %T | Ratio | 62 CPU | %T | Ratio | 4 GPU (est.) | %T | Category |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| **KSPSolve** | 521.0 | 91% | 1.0 | 114.2 | 67% | 1.0 | 53.4 (actual) | 84% | Solver total |
| **MatMult** | 435.8 | 76% | 1.0 | 94.7 | 35% | 2.3 | ~23.4 | 40% | Compute |
| **MatMultAdd** | 67.1 | 11% | 1.1 | 20.6 | 7% | 3.9 | ~3.5 | 6% | Compute |
| **KSPGMRESOrthog** | 89.3 | 15% | 1.0 | 63.9 | 28% | 3.6 | ~24.0 | 41% | Communication |
| **VecMDot** | 57.1 | 10% | 1.0 | 60.4 | 25% | **6.7** | ~24.0 | 41% | Communication |
| **VecNorm** | 13.5 | 2% | 1.1 | 11.6 | 5% | **9.4** | ~6.4 | 11% | Communication |
| **VecNormalize** | 18.4 | 3% | 1.1 | — | — | — | ~6.4 | 11% | Communication |
| **VecAXPY** | 66.7 | 11% | 1.0 | 9.2 | 3% | 3.2 | ~0.6 | 1% | Compute |
| **VecMAXPY** | 40.1 | 7% | 1.0 | 10.8 | 3% | 3.1 | ~0.0 | 0% | Compute |
| **VecScale** | 14.4 | 2% | 1.0 | 1.9 | 1% | 2.4 | ~0.0 | 0% | Compute |
| **VecCopy** | 11.3 | 2% | 1.0 | 1.9 | 1% | 2.7 | ~0.0 | 0% | — |
| **VecScatterBegin** | 5.1 | 1% | 5.0 | 8.4 | 2% | **15.6** | ~11.7 | 20% | Communication |
| **VecScatterEnd** | 15.1 | 2% | 2.1 | 12.1 | 4% | 3.5 | ~8.2 | 14% | Communication |
| **PCApply** | 10.1 | 2% | 1.0 | 1.8 | 1% | 2.6 | ~0.0 | 0% | — |

"Ratio" = max time / min time across ranks (load imbalance indicator). Higher = worse balance.

## Key Observations

Speedups below use **Quandary Used Time** (563.9s / 119.0s / 58.6s). PETSc event times from `-log_view`.

1. **GPU compute is fast, communication dominates**: On GPU, VecScatter (halo exchange) + VecMDot/VecNorm (allreduces) account for ~86% of estimated time. Actual compute (MatMult local work, VecAXPY, etc.) is negligible.

2. **62 CPU ranks hit scaling wall**: 15.5× more ranks → only 4.7× speedup (30% parallel efficiency). VecMDot ratio of 6.7 and VecNorm ratio of 9.4 show severe load imbalance from allreduce synchronization.

3. **MatMult scales well on CPU** (4.6× for 15.5× ranks = 30% efficiency) but **GMRES orthogonalization does not** (VecMDot: 57.1s → 60.4s, barely faster with 15× more ranks).

4. **GPU wins overall**: 4 GPUs beat 62 CPUs by 2.0× on Quandary Used Time despite communication overhead, because the local compute speedup is so large (~18.6× on MatMult alone).

## GPU Event Times with `-log_view_gpu_time`

A separate run with `-log_view_gpu_time` enabled gives actual per-event GPU times instead of `n/a`.
This flag forces GPU synchronization after every PETSc event, which inflates the total time
(Quandary Used Time: 86.2s vs 58.6s without the flag). The absolute times below are therefore
slower than real performance, but the **relative proportions between GPU events** and the
**CPU-to-GPU speedup ratios** are meaningful.

| Event | 4 CPU (sec) | 4 GPU (sec) | GPU speedup | GPU Mflop/s | Category |
|---|---:|---:|---:|---:|---|
| **KSPSolve** | 521.0 | 69.7 | 7.5× | 83,185 | Solver total |
| **MatMult** | 435.8 | 51.6 | 8.4× | 105,586 | Compute |
| **MatMultAdd** | 67.1 | 6.2 | 10.8× | 184,145 | Compute |
| **KSPGMRESOrthog** | 89.3 | 26.8 | 3.3× | 41,309 | Communication |
| **VecMDot** | 57.1 | 25.4 | 2.2× | 21,754 | Communication |
| **VecNorm** | 13.5 | 6.0 | 2.3× | 21,374 | Communication |
| **VecAXPY** | 66.7 | 5.3 | 12.6× | 225,183 | Compute |
| **VecMAXPY** | 40.1 | 1.7 | 23.6× | 414,658 | Compute |
| **VecScale** | 14.4 | 1.9 | 7.6× | 109,476 | Compute |
| **VecScatterBegin** | 5.1 | 9.3 | 0.5× (slower) | — | Communication |
| **VecScatterEnd** | 15.1 | 19.1 | 0.8× (slower) | — | Communication |
| **PCApply** | 10.1 | 0.6 | 16.8× | — | — |

Key findings:
- **Compute operations are 8–24× faster** on GPU (MatMult 8.4×, VecAXPY 12.6×, VecMAXPY 23.6×).
- **Allreduces (VecMDot, VecNorm) only 2.2× faster** — still latency-bound even on GPU.
- **VecScatter (halo exchange) is slower on GPU** — the GPU→MPI→GPU path adds synchronization overhead.
- **Communication is ~65% of GPU KSPSolve time**: VecMDot (25.4s) + VecScatter (28.3s) = 53.7s out of 69.7s.

Source: `benchmarks/benchmark_results_20260414_091829/` (4 CPU ranks, 4 GPU ranks, `-log_view_gpu_time` enabled).

## Configuration

- System: `nlevels = [31, 50]`, Lindblad (open system with decay)
- State dimension: 31 × 50 = 1,550; density matrix = 2,402,500 elements
- ntime = 1207, dt = 8.28e-05
- GMRES with IMR timestepper, `usematfree = false`
- Platform: Tuolumne (LLNL), AMD MI300A APUs, HPE Slingshot 11
