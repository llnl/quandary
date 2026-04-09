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

## PETSc Logging: What You Get

The runs use `-log_view -log_summary`, which provides:
- `Time (sec)` in the PETSc “Performance Summary”: PETSc total wall time **including** initialization/setup (not just the solve loop).
- Event breakdown tables including (when present) `MatMult` and `VecScatterBegin` timing percentages.

Our result extraction (`util/extract_results.sh`) reports both:
- **Quandary time**: from `timing.dat` (solve loop only; excludes setup).
- **PETSc time**: from `Time (sec):` (total runtime including setup).

Note: on GPUs, some per-event GPU timing fields can still show as `n/a` in `-log_view` depending on the backend/timer support; the top-level PETSc `Time (sec)` should still be present on successful runs.
