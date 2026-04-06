#!/usr/bin/env bash
# Benchmark Quandary CPU vs GPU performance on LC machines
set -uo pipefail

usage() {
  cat <<'EOF'
Benchmark Quandary CPU vs GPU performance using a single GPU-capable build.

Usage:
  bash util/benchmark_gpu_cpu.sh [options]

Options:
  --machine {tioga|tuolumne}  LC machine (default: tuolumne)
  --variants {cpu|gpu|both}   Which to test (default: both)
  --nprocs N                  MPI ranks (default: 8)
  --toml PATH                 Config file (default: tests/performance/configs/nlevels_32_32_32_32.toml)
  --output-dir PATH           Output directory (default: benchmark_results_<timestamp>)
  --build-only                Only build Spack environment, don't run
  --no-build                  Skip build, use existing environment (must be built first)
  -h, --help                  Show help

Note:
  Uses a single Spack environment with PETSc+Kokkos+ROCm.
  CPU mode: runs without GPU flags (uses default PETSc vectors)
  GPU mode: runs with -vec_type kokkos -mat_type aijkokkos

Examples:
  # Quick CPU vs GPU comparison (default: tuolumne, both variants, 8 ranks)
  ./util/benchmark_gpu_cpu.sh

  # Build environment only (e.g., in a batch job)
  ./util/benchmark_gpu_cpu.sh --build-only

  # Run using pre-built environment
  ./util/benchmark_gpu_cpu.sh --no-build --variants gpu --nprocs 4

  # Custom config
  ./util/benchmark_gpu_cpu.sh --toml tests/performance/configs/nlevels_16_16_16_16.toml
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# Defaults
MACHINE="tuolumne"
VARIANTS="both"
NPROCS=8
TOML="tests/performance/configs/nlevels_32_32_32_32.toml"
OUTPUT_DIR=""
BUILD_ONLY=0
NO_BUILD=0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --machine) MACHINE="$2"; shift 2;;
    --variants) VARIANTS="$2"; shift 2;;
    --nprocs) NPROCS="$2"; shift 2;;
    --toml) TOML="$2"; shift 2;;
    --output-dir) OUTPUT_DIR="$2"; shift 2;;
    --build-only) BUILD_ONLY=1; shift;;
    --no-build) NO_BUILD=1; shift;;
    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

# Validate exclusive options
if [ "$BUILD_ONLY" -eq 1 ] && [ "$NO_BUILD" -eq 1 ]; then
  die "Cannot specify both --build-only and --no-build"
fi

# Validate
command -v spack >/dev/null 2>&1 || die "spack not found in PATH"
[[ "$MACHINE" =~ ^(tioga|tuolumne)$ ]] || die "Invalid machine: $MACHINE (must be tioga or tuolumne)"
[[ "$VARIANTS" =~ ^(cpu|gpu|both)$ ]] || die "Invalid variants: $VARIANTS (must be cpu, gpu, or both)"
[[ -f "$TOML" ]] || die "Config file not found: $TOML"

# Set machine-specific defaults
case "$MACHINE" in
  tioga)
    AMDGPU_TARGET="gfx90a"
    RADIUSS_CONFIG="toss_4_x86_64_ib_cray/tioga"
    ;;
  tuolumne)
    AMDGPU_TARGET="gfx942"
    RADIUSS_CONFIG="toss_4_x86_64_ib_cray/tuolumne"
    ;;
esac

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
echo "Machine:  $MACHINE"
if [ "$BUILD_ONLY" -eq 0 ]; then
  echo "Variants: ${RUN_VARIANTS[*]}"
  echo "Ranks:    $NPROCS"
  echo "Config:   $TOML"
  echo "Output:   $OUTPUT_DIR"
fi
echo ""

# Build GPU-capable environment
if [ "$NO_BUILD" -eq 0 ]; then
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
    echo "Adding quandary@main ^petsc@3.24+kokkos+rocm amdgpu_target=${AMDGPU_TARGET} ^hipblas-common@7.2 ^hip@7.2"
    spack add "quandary@main ^petsc@3.24+kokkos+rocm amdgpu_target=${AMDGPU_TARGET} ~mmg~parmmg~saws~examples~exodusii~zoltan ^hipblas-common@7.2 ^hip@7.2"

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
    echo "  bash util/benchmark_gpu_cpu.sh --no-build"
    exit 0
  fi
  echo ""
else
  echo "Skipping build (--no-build): using existing environment"
  [ -d "$ENV_DIR" ] || die "Environment not found: $ENV_DIR (remove --no-build or build first)"
fi

# Activate environment for runs
eval $(spack env activate --sh -d "$ENV_DIR")
QUANDARY_BIN=$(spack find --format '{prefix}' quandary)/bin/quandary
[ -x "$QUANDARY_BIN" ] || die "Quandary binary not found or not executable: $QUANDARY_BIN"

echo "Binary: $QUANDARY_BIN"
echo ""

# Run each variant
for variant in "${RUN_VARIANTS[@]}"; do
  echo "========================================"
  echo "Running: $variant mode"
  echo "========================================"

  # Set PETSc runtime options
  if [ "$variant" = "cpu" ]; then
    PETSC_RUNTIME_OPTS="-log_view -log_summary"
    LAUNCHER="flux run"
    unset MPICH_GPU_SUPPORT_ENABLED
  else  # gpu
    PETSC_RUNTIME_OPTS="-vec_type kokkos -mat_type aijkokkos -use_gpu_aware_mpi 1 -log_view -log_summary"
    # Tioga has discrete GPUs (MI250X) - binding helps locality
    # Tuolumne has APUs (MI300A) - binding not supported/needed with unified memory
    if [ "$MACHINE" = "tioga" ]; then
      LAUNCHER="flux run --gpus-per-task=1 --gpu-bind=closest"
    else
      LAUNCHER="flux run --gpus-per-task=1"
    fi
    export MPICH_GPU_SUPPORT_ENABLED=1
  fi

  # Set up run directory
  RUN_DIR="${OUTPUT_DIR}/run_${variant}_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$RUN_DIR"
  cp "$TOML" "$RUN_DIR/config.toml"

  # Update config to output to run directory
  sed -i.bak "s|directory = .*|directory = \"${RUN_DIR}\"|" "$RUN_DIR/config.toml"
  rm "$RUN_DIR/config.toml.bak"

  # Run
  echo "Launcher: $LAUNCHER -n $NPROCS"
  echo "PETSc options: $PETSC_RUNTIME_OPTS"
  echo ""

  /usr/bin/time -v $LAUNCHER -n $NPROCS \
    "$QUANDARY_BIN" "$RUN_DIR/config.toml" \
    --petsc-options "$PETSC_RUNTIME_OPTS" \
    --quiet \
    > "${OUTPUT_DIR}/${variant}.log" 2>&1
  EXIT_CODE=$?

  if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ $variant completed successfully"

    # Extract key metrics
    TIME=$(grep "Used Time:" "${OUTPUT_DIR}/${variant}.log" 2>/dev/null | awk '{print $3}' || echo "N/A")
    MEMORY=$(grep "Global Memory:" "${OUTPUT_DIR}/${variant}.log" 2>/dev/null | awk '{print $3 " " $4}' || echo "N/A")
    echo "  Time: ${TIME}s"
    echo "  Memory: $MEMORY"
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
