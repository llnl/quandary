# GPU Notes: Quandary + PETSc on ROCm (Tuolumne / MI300A)

## What Works (today)

**Build:** PETSc+Kokkos+ROCm on Tuolumne (gfx942 / MI300A APU) works with Quandary using PETSc runtime options.

**Runtime requirements (important):**
- Propagate to all ranks:
  - `HSA_XNACK=1` (Kokkos/HMM on MI300A)
  - `MPICH_GPU_SUPPORT_ENABLED=1` (required for stable multi-rank GPU runs we observed)
- PETSc options for GPU:
  - `-vec_type kokkos -mat_type aijkokkos -log_view -log_summary`

The wrapper script `util/benchmark_gpu_cpu.sh` now bakes these in for GPU runs.

**Code requirement:** Quandary’s PETSc `MatShell` must have its VecType set from `-vec_type` (already implemented in `src/mastereq.cpp`).

## Build Settings (Spack pins)

These are the pins used by `util/benchmark_gpu_cpu.sh`:
- `PETSc 3.24.4` with `+kokkos+rocm amdgpu_target=gfx942`
- `Kokkos 4.6.02` with `+apu+rocm amdgpu_target=gfx942` (`cxxstd=17`)
- `ROCm/HIP 6.4.3`
- `llvm-amdgpu 6.4.3`

To force a rebuild: `rm -rf envs-benchmark/gpu_capable`

## Simple Comparisons We Run

1) **Single case (CPU vs GPU at same ranks)**
- Script: `util/benchmark_gpu_cpu.sh`
- Example: `./util/benchmark_gpu_cpu.sh --variants both --nprocs 4 --toml tests/performance/configs/nlevels_32_32_32_32.toml`

2) **Small suite (Tuolumne-only)**
- Script: `util/simple_benchmark.sh`
- Test 1: problem sizes `4^4`, `16^4`, `32^4` (CPU vs GPU at 4 ranks)
- Test 2: strong scaling for `32^4` at ranks `1,2,4` (CPU vs GPU at the same ranks)
- Test 3: fixed GPUs for `32^4` (allocate 4 GPUs and vary ranks `1,4,8`)

## Example Results (Tuolumne)

Numbers below are from `benchmark_20260409_092656` / follow-on runs.

Full `util/extract_results.sh` output:

```
================================
Test 1: Problem Size Scaling
================================

Scenario: vary problem size (`nlevels=4,16,32` → `nlevels^4` DOFs), comparing CPU vs GPU at 4 MPI ranks / 4 GPUs (1 rank per GPU).

Config  DOFs     CPU_Quandary(s)  CPU_PETSc(s)  CPU_MatMult%  CPU_Scatter%  GPU_Quandary(s)  GPU_PETSc(s)  GPU_MatMult%  GPU_Scatter%  Speedup
4^4     256      2.00867202e-01   9.209e+00     2             1             4.12218339e+00   4.615e+00     78            21            0.05x
16^4    65536    7.93781317e-01   9.655e+00     7             0             2.30298332e+00   2.854e+00     78            6             0.34x
32^4    1048576  2.44638380e+01   3.429e+01     61            1             8.13381576e+00   1.112e+01     65            8             3.01x

================================
Test 2: Strong Scaling (32^4)
================================

Scenario: fixed problem size (`nlevels=32`), compare CPU vs GPU while scaling MPI ranks and GPUs together (`ranks=GPUs=1,2,4`; 1 rank per GPU).

Ranks  CPU_Quandary(s)  CPU_PETSc(s)  CPU_MatMult%  CPU_Scatter%  GPU_Quandary(s)  GPU_PETSc(s)  GPU_MatMult%  GPU_Scatter%
1      8.35580062e+01   9.287e+01     75            0             3.87626611e+00   1.633e+01     9             0
2      4.33725220e+01   5.279e+01     70            0             1.41891114e+01   1.525e+01     53            5
4      2.52906789e+01   3.480e+01     61            1             6.77822142e+00   9.326e+00     59            9

================================
Test 3: Fixed GPUs (32^4)
================================

Scenario: fixed problem size (`nlevels=32`) and fixed GPUs (`GPUs=4` on 1 node), vary MPI ranks (`ranks=1,4,8`) to see under/over-subscription effects.

Ranks  GPUs  GPU_Quandary(s)  GPU_PETSc(s)  GPU_MatMult%  GPU_Scatter%
1      4     4.27785338e+00   1.775e+01     9             0
4      4     3.52674378e+00   6.214e+00     47            5
8      4     4.59441909e+00   6.207e+00     67            7
```

Very brief interpretation:
- Test 1: GPUs hurt at `4^4`/`16^4` (overhead dominates) but help at `32^4` (~3× solve-loop speedup).
- Test 2: GPU wins at all ranks, but the 2-rank GPU point is anomalously slow vs 1 and 4.
- Test 3: best of these is 4 ranks on 4 GPUs (1 rank/GPU); oversubscribing to 8 ranks hurts.

## PETSc Logging: What You Get

The runs use `-log_view -log_summary`, which provides:
- `Time (sec)` in the PETSc “Performance Summary”: PETSc total wall time **including** initialization/setup (not just the solve loop).
- Event breakdown tables including (when present) `MatMult` and `VecScatterBegin` timing percentages.

Our result extraction (`util/extract_results.sh`) reports both:
- **Quandary time**: from `timing.dat` (solve loop only; excludes setup).
- **PETSc time**: from `Time (sec):` (total runtime including setup).

Note: on GPUs, some per-event GPU timing fields can still show as `n/a` in `-log_view` depending on the backend/timer support; the top-level PETSc `Time (sec)` should still be present on successful runs.
