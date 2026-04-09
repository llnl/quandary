#!/usr/bin/env bash
# Benchmark Quandary CPU vs GPU performance on Tuolumne (MI300A)
set -uo pipefail

usage() {
  cat <<'EOF'
Benchmark Quandary CPU vs GPU performance using a single GPU-capable build.

Usage:
  bash util/benchmark_gpu_cpu.sh [options]

Options:
  --variants {cpu|gpu|both}   Which to test (default: both)
  --nprocs N                  MPI ranks (default: 4)
  --gpus-per-node N           GPU run: allocate N GPUs on 1 node and run NPROCS tasks (allows oversubscription)
  --toml PATH                 Config file (default: tests/performance/configs/nlevels_32_32_32_32.toml)
  --output-dir PATH           Output directory (default: benchmark_results_<timestamp>)
  --quiet                     Pass Quandary quiet flag (--quiet)
  --build-only                Only build Spack environment, don't run
  -h, --help                  Show help

Note:
  Tuolumne-only (MI300A).
  Uses a single Spack environment with PETSc+Kokkos+ROCm.
  CPU mode: runs without GPU flags (uses default PETSc vectors)
  GPU mode: runs with -vec_type kokkos -mat_type aijkokkos
  To force a rebuild: rm -rf envs-benchmark/gpu_capable

Examples:
  # Quick CPU vs GPU comparison (default: both variants, 4 ranks)
  ./util/benchmark_gpu_cpu.sh

  # Build environment only (e.g., in a batch job)
  ./util/benchmark_gpu_cpu.sh --build-only

  # Custom config
  ./util/benchmark_gpu_cpu.sh --toml tests/performance/configs/nlevels_16_16_16_16.toml
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# Defaults (tuolumne-only)
VARIANTS="both"
NPROCS=4
GPUS_PER_NODE=""
TOML="tests/performance/configs/nlevels_32_32_32_32.toml"
OUTPUT_DIR=""
BUILD_ONLY=0
QUIET=0

# Build pins (edit here if needed)
PETSC_VERSION="3.24.4"
ROCM_VERSION="6.4.3"
KOKKOS_VERSION="4.6.02"
KOKKOS_CXXSTD="17"
LLVM_AMDGPU_VERSION="6.4.3"
AMDGPU_TARGET="gfx942"
RADIUSS_CONFIG="toss_4_x86_64_ib_cray/tuolumne"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --variants) VARIANTS="$2"; shift 2;;
    --nprocs) NPROCS="$2"; shift 2;;
    --gpus-per-node) GPUS_PER_NODE="$2"; shift 2;;
    --toml) TOML="$2"; shift 2;;
    --output-dir) OUTPUT_DIR="$2"; shift 2;;
    --quiet) QUIET=1; shift;;
    --build-only) BUILD_ONLY=1; shift;;
    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

# Validate
command -v spack >/dev/null 2>&1 || die "spack not found in PATH"
[[ "$VARIANTS" =~ ^(cpu|gpu|both)$ ]] || die "Invalid variants: $VARIANTS (must be cpu, gpu, or both)"
[[ -f "$TOML" ]] || die "Config file not found: $TOML"
[[ "$KOKKOS_CXXSTD" =~ ^(14|17|20|23)$ ]] || die "Invalid KOKKOS_CXXSTD (edit near top): $KOKKOS_CXXSTD"

# Set output directory (only if we're going to run)
if [ "$BUILD_ONLY" -eq 0 ]; then
  if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="benchmark_results_$(date +%Y%m%d_%H%M%S)"
  fi
  mkdir -p "$OUTPUT_DIR"
fi

# Set which variants to run
case "$VARIANTS" in
  cpu) RUN_VARIANTS=("cpu");;
  gpu) RUN_VARIANTS=("gpu");;
  both) RUN_VARIANTS=("cpu" "gpu");;
esac

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

ENV_ROOT="${REPO_ROOT}/envs-benchmark"
ENV_DIR="${ENV_ROOT}/gpu_capable"
RADIUSS_CONFIGS="${REPO_ROOT}/.ci-scripts/radiuss-spack-configs"

[[ -d "$RADIUSS_CONFIGS" ]] || die "radiuss-spack-configs not found: $RADIUSS_CONFIGS"

echo "=== Quandary GPU/CPU Benchmark ==="
echo "Machine:  tuolumne"
if [ "$BUILD_ONLY" -eq 0 ]; then
  echo "Variants: ${RUN_VARIANTS[*]}"
  echo "Ranks:    $NPROCS"
  echo "Config:   $TOML"
  echo "Output:   $OUTPUT_DIR"
fi
echo ""

# Build GPU-capable environment
echo "========================================"
echo "Setting up GPU-capable Spack environment"
echo "========================================"

if [ -d "$ENV_DIR" ]; then
  echo "Environment exists, activating..."
  eval $(spack env activate --sh -d "$ENV_DIR")
else
  echo "Creating new environment..."
  spack env create -d "$ENV_DIR"
  eval $(spack env activate --sh -d "$ENV_DIR")

  # Add radiuss machine configs (include config.yaml and packages.yaml directly)
  spack config add "config:install_tree:padded_length:128"
  spack config add "include:[${RADIUSS_CONFIGS}/${RADIUSS_CONFIG}/config.yaml]"
  spack config add "include:[${RADIUSS_CONFIGS}/${RADIUSS_CONFIG}/packages.yaml]"

  # Add Quandary with GPU-capable PETSc
  # Tuolumne MI300A APU: enable Kokkos APU variant
  COMPILER_SPEC="%llvm-amdgpu@=${LLVM_AMDGPU_VERSION}"
  echo "Adding quandary@main (PETSc+Kokkos+ROCm, gfx942)..."
  # Note: kokkos-kernels@4.6.02 does not have +rocm/amdgpu_target variants in Spack.
  spack add "quandary@main ${COMPILER_SPEC} ^cray-mpich${COMPILER_SPEC} ^petsc@${PETSC_VERSION}+kokkos+rocm amdgpu_target=${AMDGPU_TARGET} ~mmg~parmmg~saws~examples~ml~exodusii~zoltan ${COMPILER_SPEC} ^kokkos@${KOKKOS_VERSION}+apu+rocm cxxstd=${KOKKOS_CXXSTD} amdgpu_target=${AMDGPU_TARGET} ${COMPILER_SPEC} ^kokkos-kernels@${KOKKOS_VERSION} ${COMPILER_SPEC} ^hip@${ROCM_VERSION}${COMPILER_SPEC} ^hipblas@${ROCM_VERSION}${COMPILER_SPEC} ^hipblas-common@${ROCM_VERSION}"

  # Develop mode
  spack develop -p "$REPO_ROOT" quandary@main

  echo "Concretizing..."
  spack concretize -f

  echo "Installing..."
  spack install --verbose
fi

# Find quandary binary
QUANDARY_BIN=$(spack find --format '{prefix}' quandary)/bin/quandary
[ -x "$QUANDARY_BIN" ] || die "Quandary binary not found or not executable: $QUANDARY_BIN"

echo "Binary: $QUANDARY_BIN"
eval $(spack env deactivate --sh)

if [ "$BUILD_ONLY" -eq 1 ]; then
  echo ""
  echo "=== Build Complete ==="
  echo "Binary: $QUANDARY_BIN"
  echo ""
  echo "Run benchmarks with:"
  echo "  bash util/benchmark_gpu_cpu.sh"
  exit 0
fi
echo ""

# Activate environment for runs
eval $(spack env activate --sh -d "$ENV_DIR")
QUANDARY_BIN=$(spack find --format '{prefix}' quandary)/bin/quandary
[ -x "$QUANDARY_BIN" ] || die "Quandary binary not found or not executable: $QUANDARY_BIN"

echo "Binary: $QUANDARY_BIN"
echo ""

# Optional Quandary args.
QUANDARY_ARGS=()
if [ "$QUIET" -eq 1 ]; then
  QUANDARY_ARGS+=(--quiet)
fi

# Run each variant
for variant in "${RUN_VARIANTS[@]}"; do
  echo "========================================"
  echo "Running: $variant mode"
  echo "========================================"

  # Set PETSc runtime options
  if [ "$variant" = "cpu" ]; then
    PETSC_OPTS="-log_view -log_summary"
    LAUNCHER="flux run"
    LAUNCHER_HAS_TASKS=0
  else  # gpu
    # MI300A APU: propagate XNACK + MPICH GPU support to all ranks.
    if [ -n "$GPUS_PER_NODE" ]; then
      LAUNCHER="flux run --env=HSA_XNACK=1 --env=MPICH_GPU_SUPPORT_ENABLED=1 -N 1 --gpus-per-node=${GPUS_PER_NODE} --tasks-per-node=${NPROCS}"
      LAUNCHER_HAS_TASKS=1
    else
      LAUNCHER="flux run --env=HSA_XNACK=1 --env=MPICH_GPU_SUPPORT_ENABLED=1 --gpus-per-task=1"
      LAUNCHER_HAS_TASKS=0
    fi
    PETSC_OPTS="-vec_type kokkos -mat_type aijkokkos -log_view -log_summary"
  fi

  # Set up run directory
  RUN_DIR="${OUTPUT_DIR}/run_${variant}_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$RUN_DIR"
  cp "$TOML" "$RUN_DIR/config.toml"

  # Update config to output to run directory
  sed -i.bak "s|directory = .*|directory = \"${RUN_DIR}\"|" "$RUN_DIR/config.toml"
  rm "$RUN_DIR/config.toml.bak"

  # Run
  if [ "$LAUNCHER_HAS_TASKS" -eq 1 ]; then
    echo "Command: $LAUNCHER $QUANDARY_BIN ${QUANDARY_ARGS[*]:-} $RUN_DIR/config.toml --petsc-options \"${PETSC_OPTS}\""
  else
    echo "Command: $LAUNCHER -n $NPROCS $QUANDARY_BIN ${QUANDARY_ARGS[*]:-} $RUN_DIR/config.toml --petsc-options \"${PETSC_OPTS}\""
  fi
  echo ""

  # Build command with proper quoting - use array to avoid eval quoting issues
  LAUNCHER_ARRAY=($LAUNCHER)
  if [ "$LAUNCHER_HAS_TASKS" -eq 1 ]; then
    /usr/bin/time -v "${LAUNCHER_ARRAY[@]}" \
      "$QUANDARY_BIN" "${QUANDARY_ARGS[@]}" "$RUN_DIR/config.toml" \
      --petsc-options "$PETSC_OPTS" \
      > "${OUTPUT_DIR}/${variant}.log" 2>&1
  else
    /usr/bin/time -v "${LAUNCHER_ARRAY[@]}" -n "$NPROCS" \
    "$QUANDARY_BIN" "${QUANDARY_ARGS[@]}" "$RUN_DIR/config.toml" \
    --petsc-options "$PETSC_OPTS" \
    > "${OUTPUT_DIR}/${variant}.log" 2>&1
  fi
  EXIT_CODE=$?

  if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ $variant completed successfully"

    # Extract Quandary timer from timing.dat in the run directory
    TIMING_FILE="${RUN_DIR}/timing.dat"
    if [ -f "$TIMING_FILE" ]; then
      Q_TIME=$(awk '{print $2}' "$TIMING_FILE")
      echo "  Quandary time: ${Q_TIME}s"
    fi
  else
    echo "✗ $variant failed (exit code: $EXIT_CODE)"
    echo "  Check log: ${OUTPUT_DIR}/${variant}.log"
  fi

  echo ""
done

eval $(spack env deactivate --sh)

echo "=== Benchmark Complete ==="
echo "Results in: $OUTPUT_DIR"
echo ""
echo "Extract results with:"
echo "  bash util/extract_results.sh $OUTPUT_DIR"

exit 0
