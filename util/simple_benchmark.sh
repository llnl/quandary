#!/usr/bin/env bash
# Simple CPU vs GPU benchmark suite
set -euo pipefail

MACHINE="tuolumne"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="benchmark_${TIMESTAMP}"
mkdir -p "$OUTPUT_DIR"

echo "=== Quandary CPU vs GPU Benchmark ==="
echo "Machine: $MACHINE"
echo "Output: $OUTPUT_DIR"
echo ""

# Test 1: Problem size (8 ranks, CPU vs GPU)
echo "Test 1: Problem Size Scaling (8 ranks)"
echo "======================================="
for nlevels in 4 16 32; do
  toml="tests/performance/configs/nlevels_${nlevels}_${nlevels}_${nlevels}_${nlevels}.toml"

  echo "Size: ${nlevels}^4"
  echo "  CPU..."
  ./util/benchmark_gpu_cpu.sh --machine "$MACHINE" --variants cpu \
    --nprocs 8 --toml "$toml" \
    --output-dir "${OUTPUT_DIR}/size_${nlevels}_cpu" \
    2>&1 | tee -a "${OUTPUT_DIR}/test1.log"

  echo "  GPU..."
  ./util/benchmark_gpu_cpu.sh --machine "$MACHINE" --variants gpu \
    --nprocs 8 --toml "$toml" \
    --output-dir "${OUTPUT_DIR}/size_${nlevels}_gpu" \
    2>&1 | tee -a "${OUTPUT_DIR}/test1.log"
done

echo ""
echo "Test 2: Strong Scaling (32^4 system)"
echo "====================================="
for nranks in 1 4 8; do
  toml="tests/performance/configs/nlevels_32_32_32_32.toml"

  echo "Ranks: $nranks"
  echo "  CPU..."
  ./util/benchmark_gpu_cpu.sh --machine "$MACHINE" --variants cpu \
    --nprocs $nranks --toml "$toml" \
    --output-dir "${OUTPUT_DIR}/scaling_${nranks}ranks_cpu" \
    2>&1 | tee -a "${OUTPUT_DIR}/test2.log"

  echo "  GPU..."
  ./util/benchmark_gpu_cpu.sh --machine "$MACHINE" --variants gpu \
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
