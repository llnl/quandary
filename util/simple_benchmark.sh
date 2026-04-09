#!/usr/bin/env bash
# Simple CPU vs GPU benchmark suite (tuolumne-only)
set -uo pipefail

# Optional: run only a specific test.
ONLY_TEST="all"
OUTPUT_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY_TEST="$2"; shift 2;;
    --output-dir) OUTPUT_DIR="$2"; shift 2;;
    -h|--help)
      cat <<'EOF'
Usage:
  bash util/simple_benchmark.sh [--only {1|2|3|all}] [--output-dir PATH]

Options:
  --only         Run only one test section (default: all)
  --output-dir   Reuse an existing output directory; if omitted, create benchmark_<timestamp>
EOF
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 1;;
  esac
done

# Tuolumne-only.
HOSTNAME=$(hostname)
if [[ ! "$HOSTNAME" =~ tuolumne ]]; then
  echo "ERROR: util/simple_benchmark.sh is tuolumne-only (hostname: ${HOSTNAME})" >&2
  exit 1
fi
MACHINE="tuolumne"
MAX_GPUS=4

if [ -z "$OUTPUT_DIR" ]; then
  TIMESTAMP=$(date +%Y%m%d_%H%M%S)
  OUTPUT_DIR="benchmark_${TIMESTAMP}"
fi
mkdir -p "$OUTPUT_DIR"

echo "=== Quandary CPU vs GPU Benchmark ==="
echo "Machine: $MACHINE"
echo "Max GPUs: $MAX_GPUS"
echo "Output: $OUTPUT_DIR"
echo ""

# Test 1: Problem size (compare CPU vs GPU at same ranks)
if [ "$ONLY_TEST" = "all" ] || [ "$ONLY_TEST" = "1" ]; then
  echo "Test 1: Problem Size Scaling"
  echo "======================================="
  for nlevels in 4 16 32; do
    toml="tests/performance/configs/nlevels_${nlevels}_${nlevels}_${nlevels}_${nlevels}.toml"

    echo "Size: ${nlevels}^4"
    echo "  CPU ($MAX_GPUS ranks)..."
    ./util/benchmark_gpu_cpu.sh --variants cpu \
      --nprocs $MAX_GPUS --toml "$toml" \
      --output-dir "${OUTPUT_DIR}/size_${nlevels}_cpu" \
      2>&1 | tee -a "${OUTPUT_DIR}/test1.log"
    echo "  [DEBUG] CPU test completed with exit code: $?"

    echo "  GPU ($MAX_GPUS ranks)..."
    ./util/benchmark_gpu_cpu.sh --variants gpu \
      --nprocs $MAX_GPUS --toml "$toml" \
      --output-dir "${OUTPUT_DIR}/size_${nlevels}_gpu" \
      2>&1 | tee -a "${OUTPUT_DIR}/test1.log"
    echo "  [DEBUG] GPU test completed with exit code: $?"
  done
fi

echo ""
echo "Test 2: Strong Scaling (32^4 system, up to $MAX_GPUS GPUs)"
echo "====================================="
# Scale from 1 up to MAX_GPUS (powers of 2)
SCALING_RANKS=(1)
if [ $MAX_GPUS -ge 2 ]; then SCALING_RANKS+=(2); fi
if [ $MAX_GPUS -ge 4 ]; then SCALING_RANKS+=(4); fi
if [ $MAX_GPUS -ge 8 ]; then SCALING_RANKS+=(8); fi

if [ "$ONLY_TEST" = "all" ] || [ "$ONLY_TEST" = "2" ]; then
  for nranks in "${SCALING_RANKS[@]}"; do
    toml="tests/performance/configs/nlevels_32_32_32_32.toml"

    echo "Ranks: $nranks"
    echo "  CPU..."
    ./util/benchmark_gpu_cpu.sh --variants cpu \
      --nprocs $nranks --toml "$toml" \
      --output-dir "${OUTPUT_DIR}/scaling_${nranks}ranks_cpu" \
      2>&1 | tee -a "${OUTPUT_DIR}/test2.log"

    echo "  GPU..."
    ./util/benchmark_gpu_cpu.sh --variants gpu \
      --nprocs $nranks --toml "$toml" \
      --output-dir "${OUTPUT_DIR}/scaling_${nranks}ranks_gpu" \
      2>&1 | tee -a "${OUTPUT_DIR}/test2.log"
  done
fi

echo ""
echo "Test 3: Fixed GPUs (32^4, allocate 4 GPUs and vary ranks)"
echo "====================================="
if [ "$ONLY_TEST" = "all" ] || [ "$ONLY_TEST" = "3" ]; then
  toml="tests/performance/configs/nlevels_32_32_32_32.toml"
  FIXED_GPU_RANKS=(1 4 8)
  for nranks in "${FIXED_GPU_RANKS[@]}"; do
    echo "Ranks: $nranks (GPUs fixed at $MAX_GPUS)"
    ./util/benchmark_gpu_cpu.sh --variants gpu \
      --gpus-per-node $MAX_GPUS \
      --nprocs $nranks --toml "$toml" \
      --output-dir "${OUTPUT_DIR}/fixed_${MAX_GPUS}gpus_${nranks}ranks_gpu" \
      2>&1 | tee -a "${OUTPUT_DIR}/test3.log"
  done
fi

echo ""
echo "=== Benchmark Complete ==="
echo "Results saved to: $OUTPUT_DIR/"
echo ""
echo "Extract results with:"
echo "  bash util/extract_results.sh $OUTPUT_DIR"
