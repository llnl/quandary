#!/usr/bin/env bash
# Benchmark Quandary CPU vs GPU performance on Tuolumne (MI300A)
set -uo pipefail

usage() {
  cat <<'EOF'
Benchmark Quandary CPU vs GPU performance using a single GPU-capable build.

Usage:
  bash util/benchmark_gpu_cpu.sh <config.toml> [options]

Arguments:
  config.toml                 Quandary config file (required)

Options:
  --variants {cpu|gpu|both}   Which to test (default: both)
  --nprocs N                  MPI ranks (default: 4)
  --output-dir PATH           Output directory (default: benchmarks/benchmark_results_<timestamp>)
  --quiet                     Suppress Quandary output and PETSc logging
  -h, --help                  Show help

Note:
  Tuolumne-only (MI300A).
  Uses .spack_env_tuolumne/ (PETSc+Kokkos+ROCm).
  CPU mode: runs without GPU flags (uses default PETSc vectors)
  GPU mode: runs with -vec_type kokkos -mat_type aijkokkos

Examples:
  # CPU vs GPU comparison
  ./util/benchmark_gpu_cpu.sh tests/performance/configs/nlevels_32_32_32_32.toml

  # GPU only
  ./util/benchmark_gpu_cpu.sh tests/performance/configs/nlevels_16_16_16_16.toml --variants gpu
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# Defaults
VARIANTS="both"
NPROCS=4
TOML=""
OUTPUT_DIR=""
QUIET=0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --variants) VARIANTS="$2"; shift 2;;
    --nprocs) NPROCS="$2"; shift 2;;
    --output-dir) OUTPUT_DIR="$2"; shift 2;;
    --quiet) QUIET=1; shift;;
    -h|--help) usage; exit 0;;
    -*) die "Unknown option: $1";;
    *) TOML="$1"; shift;;
  esac
done

# Validate
[[ -n "$TOML" ]] || die "Config file required. Usage: bash util/benchmark_gpu_cpu.sh <config.toml>"
[[ -f "$TOML" ]] || die "Config file not found: $TOML"
command -v spack >/dev/null 2>&1 || die "spack not found in PATH"
[[ "$VARIANTS" =~ ^(cpu|gpu|both)$ ]] || die "Invalid variants: $VARIANTS (must be cpu, gpu, or both)"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# Set output directory
if [ -z "$OUTPUT_DIR" ]; then
  OUTPUT_DIR="benchmarks/benchmark_results_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUTPUT_DIR"

# Set which variants to run
case "$VARIANTS" in
  cpu) RUN_VARIANTS=("cpu");;
  gpu) RUN_VARIANTS=("gpu");;
  both) RUN_VARIANTS=("cpu" "gpu");;
esac

echo "=== Quandary GPU/CPU Benchmark ==="
echo "Variants: ${RUN_VARIANTS[*]}"
echo "Ranks:    $NPROCS"
echo "Config:   $TOML"
echo "Output:   $OUTPUT_DIR"
echo ""

# Activate environment and build if needed
echo "Setting up Spack environment..."
eval $(spack env activate --sh -d "${REPO_ROOT}/.spack_env_tuolumne")
spack install

QUANDARY_BIN=$(spack find --format '{prefix}' quandary)/bin/quandary
[ -x "$QUANDARY_BIN" ] || die "Quandary binary not found: $QUANDARY_BIN"
echo "Binary: $QUANDARY_BIN"
echo ""

# Run each variant
for variant in "${RUN_VARIANTS[@]}"; do
  echo "========================================"
  echo "Running: $variant mode"
  echo "========================================"

  # Set PETSc runtime options and launcher
  if [ "$variant" = "cpu" ]; then
    PETSC_OPTS=""
    LAUNCHER=(flux run -n "$NPROCS")
  else
    PETSC_OPTS="-vec_type kokkos -mat_type aijkokkos"
    LAUNCHER=(flux run --env=HSA_XNACK=1 --env=MPICH_GPU_SUPPORT_ENABLED=1 --gpus-per-task=1 -n "$NPROCS")
  fi

  # Add PETSc logging unless quiet
  if [ "$QUIET" -eq 0 ]; then
    PETSC_OPTS="$PETSC_OPTS -log_view -log_summary"
  fi

  # Quandary args
  QUANDARY_ARGS=()
  if [ "$QUIET" -eq 1 ]; then
    QUANDARY_ARGS+=(--quiet)
  fi

  # Set up run directory
  RUN_DIR="${OUTPUT_DIR}/run_${variant}_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$RUN_DIR"
  cp "$TOML" "$RUN_DIR/config.toml"
  sed -i.bak "s|directory = .*|directory = \"${RUN_DIR}\"|" "$RUN_DIR/config.toml"
  rm "$RUN_DIR/config.toml.bak"

  echo "Command: ${LAUNCHER[*]} $QUANDARY_BIN ${QUANDARY_ARGS[*]:-} $RUN_DIR/config.toml --petsc-options \"${PETSC_OPTS}\""
  echo ""

  # Run
  "${LAUNCHER[@]}" \
    "$QUANDARY_BIN" "${QUANDARY_ARGS[@]}" "$RUN_DIR/config.toml" \
    --petsc-options "$PETSC_OPTS" \
    > "${OUTPUT_DIR}/${variant}.log" 2>&1
  EXIT_CODE=$?

  if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ $variant completed successfully"
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
