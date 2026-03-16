# GPU Notes: Quandary + PETSc on ROCm (Tioga/Tuolumne)

## Summary

**Result:** Quandary achieves **11.3x speedup** on AMD MI250X GPUs (Tioga) using PETSc 3.24.4 with Kokkos backend and GPU-aware MPI.

**Key findings:**
- PETSc+Kokkos backend (`-vec_type kokkos -mat_type aijkokkos`) works well on AMD GPUs
- GPU-aware MPI is critical: provides 3.6x additional speedup over GPU-only acceleration
  - Without it: ~17,800 host↔device transfers per run (~8 GB each direction), high memory bandwidth usage
- Required fix: MatShell VecType configuration (already in `src/mastereq.cpp`)

## Performance Results (Tioga, 8 MPI ranks, 32⁴ system)

| Configuration | Time (s) | Speedup vs CPU | Memory (GB) |
|---------------|----------|----------------|-------------|
| CPU baseline | 44.3 | 1.0x | 2.7 (~340 MB/rank) |
| Kokkos (no GPU-aware MPI) | 14.1 | 3.1x | 20.0 (~2.5 GB/rank) |
| Kokkos (GPU-aware MPI) | 3.9 | **11.3x** | 20.0 (~2.5 GB/rank) |

**Where time is spent (from `-log_view`):**
- CPU: MatMult ≈ 85% of runtime
- Kokkos (no GPU-aware): MatMult ≈ 83% (communication-bound)
- Kokkos (GPU-aware): MatMult ≈ 52% (compute-bound)

GPU-aware MPI eliminates host staging (`GPU → host → MPI → host → GPU` becomes `GPU → MPI → GPU`), reducing communication overhead and shifting the bottleneck back to compute.

**Memory usage:** GPU runs require ~7.4x more memory due to device allocations, Kokkos dual-view structures, and PETSc matrix format conversions during setup.

**PETSc configuration:**
- CPU: `^petsc@3.24.4~kokkos~rocm`
- GPU: `^petsc@3.24.4+kokkos+rocm amdgpu_target=gfx90a`
- Runtime: `-vec_type kokkos -mat_type aijkokkos -use_gpu_aware_mpi 1`

## Background

- Quandary supports passing PETSc runtime flags via `quandary <config> --petsc-options “<...>”` (no user code changes needed).
- The main operator is a PETSc `MatShell` (`src/mastereq.cpp`); time stepping repeatedly calls `MatMult` and `KSPSolve` (`src/timestepper.cpp`).
- **PETSc backend support on AMD GPUs:**
  - **Kokkos (recommended):** Vec and Mat fully supported, production-ready
  - **HIP-only (`rocm`):** Vec supported, Mat in development (not recommended for production)
- PETSc GPU roadmap: https://petsc.org/release/overview/gpu_roadmap/
- **Profiling note:** PETSc's `-log_view` timers are an effective way to study Quandary performance, as ~91% of solve time is spent in PETSc operations (KSPSolve, MatMult, vector operations).

**Why Kokkos over HIP-only?** Per PETSc 3.24.4 docs, the HIP backend's matrix operations are still experimental. We tested the HIP-only path (`+rocm~kokkos`, `-vec_type hip -mat_type aijhipsparse`) and observed significant **slowdowns vs CPU** (260-326s vs 44s baseline), confirming the experimental Mat support is not production-ready. The Kokkos path provides full Vec+Mat GPU support and delivered the 11.3x speedup shown above.

## Quick Start

Use `util/gpu_petsc_ab.sh` to reproduce the benchmark or test your own configurations:

```bash
# Tioga (default: cpu,kokkos,rocm variants, 8 MPI ranks)
./util/gpu_petsc_ab.sh --machine tioga

# Tuolumne
./util/gpu_petsc_ab.sh --machine tuolumne

# Only CPU + Kokkos comparison
./util/gpu_petsc_ab.sh --variants cpu,kokkos
```

Each run creates a unique timestamped directory under `envs-gpu-ab/runs/` with logs and PETSc `-log_view` output.

## Script Details

The `util/gpu_petsc_ab.sh` script:
- Creates separate Spack environments for CPU/Kokkos/ROCm variants
- Uses RADIUSS Tioga/Tuolumne machine configs for external packages
- Defaults: PETSc 3.24.4, ROCm 6.4.3, llvm-amdgpu@6.4.3 compiler
- GPU variants: Enables GPU-aware MPI by default (`-use_gpu_aware_mpi 1` + `MPICH_GPU_SUPPORT_ENABLED=1`)
- Launcher: `flux run --gpus-per-task=1 -n 8` (Tioga) or `-n 4` (Tuolumne); Flux uses `mpibind` by default to map ranks to GPUs/CPUs
- Output: Per-run timestamped directories under `envs-gpu-ab/runs/` with logs and PETSc performance data

**PETSc options for each variant:**
- `cpu`: `-log_view -log_summary`
- `kokkos`: `-vec_type kokkos -mat_type aijkokkos -use_gpu_aware_mpi 1 -log_view -log_summary`
- `rocm`: `-vec_type hip -mat_type aijhipsparse -use_gpu_aware_mpi 1 -log_view -log_summary`
  - *Note:* HIP Mat support is experimental (in development per PETSc roadmap). Use Kokkos for production.

## Interpreting Results

**In PETSc `-log_view` output:**
- Check `KSPSolve` time (main linear solver): Should dominate for good GPU candidates
- Look for `GPU %F` column: 100% means operation ran entirely on GPU
- For Kokkos runs, confirm `GPU Mflop/s` values appear for Vec/Mat operations
- Compare `MatMult` and `VecScatterBegin/End` times: High scatter time suggests MPI overhead

**Performance indicators:**
- Good: KSPSolve dominates time, high GPU Mflop/s, GPU %F = 100%
- Poor: “User” time high, low GPU activity, or excessive VecScatter overhead

**If speedup is disappointing:**
- Check GPU-aware MPI is enabled (huge impact: 3.6x in our tests)
- Verify PETSc sees GPUs: Look for “HIP architecture 90” in header
- Profile with `rocprof` or Omniperf for detailed GPU kernel analysis

## Technical Details

### MatShell VecType Fix (Required for Kokkos)

**Problem:** PETSc `-mat_type aijkokkos` originally crashed with MatShell operators.

**Root cause:** PETSc’s AIJKokkos MatMult expects all vectors to be `VECKOKKOS`. Quandary’s MatShell operator creates work vectors using the Mat’s default VecType, which was `seq`/`mpi` even when `-vec_type kokkos` was specified. This caused type mismatches and crashes.

**Solution:** Set MatShell VecType from the runtime `-vec_type` option ([src/mastereq.cpp:102-109](src/mastereq.cpp#L102-L109)):
```cpp
char vec_type_opt[128] = {0};
PetscBool vec_type_set = PETSC_FALSE;
PetscOptionsGetString(NULL, NULL, "-vec_type", vec_type_opt, sizeof(vec_type_opt), &vec_type_set);
if (vec_type_set) MatSetVecType(RHS, vec_type_opt);
```

This ensures PETSc’s `MatCreateVecs()` creates work vectors matching the requested backend.

### Spack Tips

**CPU variant:** Don’t add ROCm-only constraints like `^hip@...` or `amdgpu_target=...` — they can implicitly force `+rocm` in dependencies and cause conflicts.

**Debug builds:** Add `--petsc-debug` to build PETSc with assertions and debug symbols:
```bash
./util/gpu_petsc_ab.sh --variants kokkos --petsc-debug --nprocs 1 \
  --debugger "gdb -batch -ex run -ex bt --args"
```

### Spack Package Improvements

**Current limitation:** The Quandary Spack package’s `+rocm` variant currently only does `^petsc+rocm`, which enables HIP-only support (Vec only, experimental Mat). Based on this testing, we should update the package:

**Proposed changes:**
- `quandary+rocm` should imply `^petsc+rocm+kokkos` (not just `+rocm`)
- Add `^kokkos` dependency when `+rocm` is enabled
- Ensure `amdgpu_target` propagates consistently to PETSc, Kokkos, and Kokkos-kernels
- Document in package.py that GPU-aware MPI runtime flags are needed for best performance

This would make `spack install quandary+rocm` work out-of-the-box with production-ready Kokkos backend instead of experimental HIP-only, and users would only need to add runtime flags (`--petsc-options "-vec_type kokkos -mat_type aijkokkos -use_gpu_aware_mpi 1"`) rather than rebuilding with custom PETSc specs.

### Comparing 64 MPI Ranks (CPU vs Kokkos)

To compare 64 MPI ranks on a single Tioga node:

**CPU (no GPUs requested):**
```bash
./util/gpu_petsc_ab.sh --variants cpu --nprocs 64 --launcher "flux run"
```

**Kokkos (64 ranks sharing 8 GPUs):**
```bash
./util/gpu_petsc_ab.sh --variants kokkos --nprocs 64 \
  --launcher "flux run --setopt=mpibind=verbose:1"
```

If you do not pass `--launcher`, `util/gpu_petsc_ab.sh` will auto-add `--gpu-bind=closest` when `flux run --help` reports that option; it’s optional and not available on all systems.

Add `--run-only` to reuse existing env installs, and `--no-quiet` to include Quandary’s own `Used Time:` line in the log.

**Observed results (Tioga, 1 node, 64 ranks; 8 GPUs shared by 64 ranks):**
- CPU log: `envs-gpu-ab/runs/cpu__20260316_085304/cpu.log`
  - `/usr/bin/time -p` `real`: ~7.27s mean across ranks
  - Quandary `Used Time:` ~5.83s
  - PETSc `Time (sec)`: ~6.02s; `MatMult`: ~5.08s
- Kokkos log: `envs-gpu-ab/runs/kokkos_vec-kokkos_mat-aijkokkos__20260316_085424/kokkos.log`
  - Flux `mpibind` maps ~8 ranks/GPU (ranks 0–7 on one GPU, 8–15 on the next, etc.)
  - `/usr/bin/time -p` `real`: ~39.17s mean across ranks (≈5.4× slower than CPU for this run)
  - Quandary `Used Time:` ~4.46s (faster than CPU *inside Quandary’s timed region*)
  - PETSc confirms the intended GPU path is enabled (`-vec_type kokkos -mat_type aijkokkos -use_gpu_aware_mpi 1`) and reports `GPU %F=100` for the main compute events, but many GPU event timings show as `n/a` in `-log_view`.

**Takeaway:** the earlier speedups in this note were measured with low MPI rank counts (e.g., 1 rank/GPU). For 64 ranks sharing 8 GPUs, the wall-clock time is much worse even though the in-code `Used Time` is lower; this suggests substantial overhead outside Quandary’s timed region (e.g., GPU runtime initialization, device setup, synchronization, or GPU oversubscription effects).
