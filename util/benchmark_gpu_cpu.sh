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
  --nprocs N                  MPI ranks (default: 8 on tioga, 4 on tuolumne)
  --toml PATH                 Config file (default: tests/performance/configs/nlevels_32_32_32_32.toml)
  --output-dir PATH           Output directory (default: benchmark_results_<timestamp>)
  --quiet                     Pass Quandary quiet flag (--quiet)
  --petsc-version VER         PETSc version (default: 3.24.4)
  --rocm-version VER          ROCm/HIP version (default: 6.4.3)
  --kokkos-version VER        Kokkos version (default: 4.6.02)
  --kokkos-cxxstd N           Kokkos C++ standard (default: 17)
  --llvm-amdgpu VER           llvm-amdgpu compiler version (default: 6.4.3)
  --build-only                Only build Spack environment, don't run
  --no-build                  Skip build, use existing environment (must be built first)
  --rebuild                   Remove existing environment and rebuild from scratch
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
NPROCS_SET=0
TOML="tests/performance/configs/nlevels_32_32_32_32.toml"
OUTPUT_DIR=""
PETSC_VERSION="3.24.4"
ROCM_VERSION="6.4.3"
KOKKOS_VERSION="4.6.02"
KOKKOS_CXXSTD="17"
LLVM_AMDGPU_VERSION="6.4.3"
BUILD_ONLY=0
NO_BUILD=0
REBUILD=0
QUIET=0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --machine) MACHINE="$2"; shift 2;;
    --variants) VARIANTS="$2"; shift 2;;
    --nprocs) NPROCS="$2"; NPROCS_SET=1; shift 2;;
    --toml) TOML="$2"; shift 2;;
    --output-dir) OUTPUT_DIR="$2"; shift 2;;
    --quiet) QUIET=1; shift;;
    --petsc-version) PETSC_VERSION="$2"; shift 2;;
    --rocm-version) ROCM_VERSION="$2"; shift 2;;
    --kokkos-version) KOKKOS_VERSION="$2"; shift 2;;
    --kokkos-cxxstd) KOKKOS_CXXSTD="$2"; shift 2;;
    --llvm-amdgpu) LLVM_AMDGPU_VERSION="$2"; shift 2;;
    --build-only) BUILD_ONLY=1; shift;;
    --no-build) NO_BUILD=1; shift;;
    --rebuild) REBUILD=1; shift;;
    -h|--help) usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

# Validate exclusive options
if [ "$BUILD_ONLY" -eq 1 ] && [ "$NO_BUILD" -eq 1 ]; then
  die "Cannot specify both --build-only and --no-build"
fi
if [ "$REBUILD" -eq 1 ] && [ "$NO_BUILD" -eq 1 ]; then
  die "Cannot specify both --rebuild and --no-build"
fi

# Validate
command -v spack >/dev/null 2>&1 || die "spack not found in PATH"
[[ "$MACHINE" =~ ^(tioga|tuolumne)$ ]] || die "Invalid machine: $MACHINE (must be tioga or tuolumne)"
[[ "$VARIANTS" =~ ^(cpu|gpu|both)$ ]] || die "Invalid variants: $VARIANTS (must be cpu, gpu, or both)"
[[ -f "$TOML" ]] || die "Config file not found: $TOML"
[[ "$KOKKOS_CXXSTD" =~ ^(14|17|20|23)$ ]] || die "Invalid --kokkos-cxxstd: $KOKKOS_CXXSTD"

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

# Machine-specific default ranks (only if user didn't pass --nprocs).
if [ "$NPROCS_SET" -eq 0 ]; then
  if [ "$MACHINE" = "tuolumne" ]; then
    NPROCS=4
  else
    NPROCS=8
  fi
fi

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
echo "PETSc:    ${PETSC_VERSION}"
echo "ROCm:     ${ROCM_VERSION}"
echo "Kokkos:   ${KOKKOS_VERSION} (cxxstd=${KOKKOS_CXXSTD})"
echo "Compiler: llvm-amdgpu@=${LLVM_AMDGPU_VERSION}"
echo ""

# Build GPU-capable environment
if [ "$NO_BUILD" -eq 0 ]; then
  echo "========================================"
  echo "Setting up GPU-capable Spack environment"
  echo "========================================"

  # Handle rebuild if requested
  if [ "$REBUILD" -eq 1 ] && [ -d "$ENV_DIR" ]; then
    echo "Removing existing environment for rebuild..."
    rm -rf "$ENV_DIR"
  fi

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
    # For MI300A APUs (tuolumne), enable Kokkos APU variant
    COMPILER_SPEC="%llvm-amdgpu@=${LLVM_AMDGPU_VERSION}"
    if [ "$MACHINE" = "tuolumne" ]; then
      echo "Adding quandary@main ^petsc@${PETSC_VERSION}+kokkos+rocm ^kokkos@${KOKKOS_VERSION}+apu+rocm cxxstd=${KOKKOS_CXXSTD} amdgpu_target=${AMDGPU_TARGET} ^hip@${ROCM_VERSION}"
      # Note: kokkos-kernels@4.6.02 does not have +rocm/amdgpu_target variants in Spack.
      spack add "quandary@main ${COMPILER_SPEC} ^cray-mpich${COMPILER_SPEC} ^petsc@${PETSC_VERSION}+kokkos+rocm amdgpu_target=${AMDGPU_TARGET} ~mmg~parmmg~saws~examples~ml~exodusii~zoltan ${COMPILER_SPEC} ^kokkos@${KOKKOS_VERSION}+apu+rocm cxxstd=${KOKKOS_CXXSTD} amdgpu_target=${AMDGPU_TARGET} ${COMPILER_SPEC} ^kokkos-kernels@${KOKKOS_VERSION} ${COMPILER_SPEC} ^hip@${ROCM_VERSION}${COMPILER_SPEC} ^hipblas@${ROCM_VERSION}${COMPILER_SPEC} ^hipblas-common@${ROCM_VERSION}"
    else
      echo "Adding quandary@main ^petsc@${PETSC_VERSION}+kokkos+rocm amdgpu_target=${AMDGPU_TARGET} ^kokkos@${KOKKOS_VERSION}+rocm cxxstd=${KOKKOS_CXXSTD} ^hip@${ROCM_VERSION}"
      spack add "quandary@main ${COMPILER_SPEC} ^cray-mpich${COMPILER_SPEC} ^petsc@${PETSC_VERSION}+kokkos+rocm amdgpu_target=${AMDGPU_TARGET} ~mmg~parmmg~saws~examples~ml~exodusii~zoltan ${COMPILER_SPEC} ^kokkos@${KOKKOS_VERSION}+rocm cxxstd=${KOKKOS_CXXSTD} amdgpu_target=${AMDGPU_TARGET} ${COMPILER_SPEC} ^kokkos-kernels@${KOKKOS_VERSION} ${COMPILER_SPEC} ^hip@${ROCM_VERSION}${COMPILER_SPEC} ^hipblas@${ROCM_VERSION}${COMPILER_SPEC} ^hipblas-common@${ROCM_VERSION}"
    fi

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
    unset MPICH_GPU_SUPPORT_ENABLED
  else  # gpu
    ENV_FLAGS=""
    # Tioga: discrete GPUs (MI250X) - use GPU-aware MPI and binding
    # Tuolumne: APUs (MI300A) - unified memory, no GPU-aware MPI needed
    if [ "$MACHINE" = "tioga" ]; then
      ENV_FLAGS="--env=MPICH_GPU_SUPPORT_ENABLED=1"
      LAUNCHER="flux run ${ENV_FLAGS} --gpus-per-task=1 --gpu-bind=closest"
      PETSC_OPTS="-vec_type kokkos -mat_type aijkokkos -use_gpu_aware_mpi 1 -log_view -log_summary"
      export MPICH_GPU_SUPPORT_ENABLED=1
    else  # tuolumne
      # MI300A APU requires XNACK for unified memory access.
      # Explicitly propagate to Flux ranks, along with GPU-aware MPI enablement.
      ENV_FLAGS="--env=HSA_XNACK=1 --env=MPICH_GPU_SUPPORT_ENABLED=1"
      LAUNCHER="flux run ${ENV_FLAGS} --gpus-per-task=1"
      PETSC_OPTS="-vec_type kokkos -mat_type aijkokkos -log_view -log_summary"
      export HSA_XNACK=1
      export MPICH_GPU_SUPPORT_ENABLED=1
    fi
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
  echo "Command: $LAUNCHER -n $NPROCS $QUANDARY_BIN ${QUANDARY_ARGS[*]:-} $RUN_DIR/config.toml --petsc-options \"${PETSC_OPTS}\""
  echo ""

  # Build command with proper quoting - use array to avoid eval quoting issues
  LAUNCHER_ARRAY=($LAUNCHER)
  /usr/bin/time -v "${LAUNCHER_ARRAY[@]}" -n "$NPROCS" \
    "$QUANDARY_BIN" "${QUANDARY_ARGS[@]}" "$RUN_DIR/config.toml" \
    --petsc-options "$PETSC_OPTS" \
    > "${OUTPUT_DIR}/${variant}.log" 2>&1
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
