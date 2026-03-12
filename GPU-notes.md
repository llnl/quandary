# GPU Notes: Quandary + PETSc on ROCm (Tioga/Tuolumne)

## Findings

- Quandary already supports passing PETSc runtime flags via `quandary <config> --petsc-options "<...>"` (no code changes needed to try PETSc GPU backends).
- The main operator is a PETSc `MatShell` (`src/mastereq.cpp`), and time stepping repeatedly calls `MatMult` and `KSPSolve` (`src/timestepper.cpp`), so GPU benefit depends on whether time is dominated by PETSc kernels vs Quandary’s shell callback / synchronization.
- PETSc docs suggest on AMD:
  - HIP/ROCm: Vec supported; Mat “in development” (may limit speedups depending on matrix path).
  - Kokkos: Vec and Mat supported (best first target on AMD).

## Recommended first experiment (simple)

- Start with wall-clock + PETSc logging (`-log_view -log_summary`); only add Caliper if PETSc logs don’t explain the result.
- Compare 2–3 PETSc builds using the same input and MPI/GPU binding:
  - `cpu`: `^petsc@3.24.4~kokkos~rocm`
  - `kokkos`: `^petsc@3.24.4+kokkos+rocm` (preferred GPU path on AMD)
  - `rocm`: `^petsc@3.24.4+rocm~kokkos` (experimental HIP-only comparison)

## Script

Use `util/gpu_petsc_ab.sh` to create Spack envs (using RADIUSS Tioga/Tuolumne configs), build, and run the A/B/C commands with PETSc logging.

Examples:

- Tioga, run everything (defaults to `flux run --gpus-per-task=1 --gpu-bind=closest` and `-n 8`):
  - `bash util/gpu_petsc_ab.sh --machine tioga`
- Tuolumne:
  - `bash util/gpu_petsc_ab.sh --machine tuolumne`
- Only CPU + Kokkos:
  - `bash util/gpu_petsc_ab.sh --variants cpu,kokkos`

The logs are written under `envs-gpu-ab/` by default (`cpu.log`, `kokkos.log`, `rocm.log`).

## Interpreting results

- In `-log_view`, check whether time is dominated by `MatMult`/`KSPSolve` (better GPU candidate) vs “user”/shell work (GPU may not help without refactors).
- If Kokkos/HIP is engaged, confirm `-vec_type`/`-mat_type` selections appear and that PETSc reports meaningful time in GPU-capable kernels.

## Next steps

1. Run `cpu` and `kokkos` first; only run `rocm` if you want the additional HIP-only comparison.
2. If GPU speedup is small and PETSc logs indicate shell/user time dominates or GPU activity is minimal, add a short Caliper run to diagnose GPU activity and sync/transfer overhead.

## Common concretization pitfall (CPU baseline)

- For the `cpu` variant, do **not** add ROCm-only constraints like `^hip@...` or `amdgpu_target=...`. Those can implicitly force `+rocm` in some dependency and create `~rocm`/`+rocm` conflicts (or `amdgpu_target=none` conflicts).
