#!/usr/bin/env bash
# Simple CPU vs GPU benchmark suite (runs on Tuolumne by default)
set -uo pipefail

# Tuolumne-only.
HOSTNAME=$(hostname)
if [[ ! "$HOSTNAME" =~ tuolumne ]]; then
  echo "ERROR: util/simple_benchmark.sh is tuolumne-only (hostname: ${HOSTNAME})" >&2
  exit 1
fi
MACHINE="tuolumne"
MAX_GPUS=4

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="benchmark_${TIMESTAMP}"
mkdir -p "$OUTPUT_DIR"

echo "=== Quandary CPU vs GPU Benchmark ==="
echo "Machine: $MACHINE"
echo "Max GPUs: $MAX_GPUS"
echo "Output: $OUTPUT_DIR"
echo ""

# Test 1: Problem size (compare CPU vs GPU at same ranks)
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

echo ""
echo "Test 2: Strong Scaling (32^4 system, up to $MAX_GPUS GPUs)"
echo "====================================="
# Scale from 1 up to MAX_GPUS (powers of 2)
SCALING_RANKS=(1)
if [ $MAX_GPUS -ge 2 ]; then SCALING_RANKS+=(2); fi
if [ $MAX_GPUS -ge 4 ]; then SCALING_RANKS+=(4); fi
if [ $MAX_GPUS -ge 8 ]; then SCALING_RANKS+=(8); fi

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

echo ""
echo "=== Benchmark Complete ==="
echo "Results saved to: $OUTPUT_DIR/"
echo ""
echo "Extract results with:"
echo "  bash util/extract_results.sh $OUTPUT_DIR"
