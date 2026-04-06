#!/usr/bin/env bash
# Extract key metrics from benchmark logs
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <benchmark_dir>"
  exit 1
fi

BENCH_DIR="$1"

extract_from_run() {
  local bench_dir="$1"
  local variant="$2"  # cpu or gpu

  local log="${bench_dir}/${variant}.log"

  if [ ! -f "$log" ]; then
    echo "n/a,n/a,n/a,n/a,n/a"
    return
  fi

  # Find the run directory (run_cpu_* or run_gpu_*)
  local run_subdir=$(find "$bench_dir" -maxdepth 1 -type d -name "run_${variant}_*" | head -n 1)

  # Quandary timer from timing.dat
  local q_time="n/a"
  if [ -n "$run_subdir" ] && [ -f "$run_subdir/timing.dat" ]; then
    q_time=$(awk '{print $2}' "$run_subdir/timing.dat")
  fi

  # DOFs from config_log.toml
  local dofs="n/a"
  if [ -n "$run_subdir" ] && [ -f "$run_subdir/config_log.toml" ]; then
    # Extract nlevels array and calculate product
    local nlevels=$(grep "^nlevels = " "$run_subdir/config_log.toml" | sed 's/.*\[\(.*\)\].*/\1/' | tr -d ' ')
    if [ -n "$nlevels" ]; then
      dofs=$(python3 -c "import math; print(math.prod([int(x) for x in '$nlevels'.split(',')]))" 2>/dev/null || echo "n/a")
    fi
  fi

  # PETSc total time (in the summary section, column 3 is Avg time)
  local p_time=$(grep "^Time (sec):" "$log" | awk '{print $3}')

  # MatMult %T (column 11 is %T in Global section)
  local matmult=$(grep "^MatMult " "$log" | awk '{print $11}')

  # VecScatterBegin %T (column 11 is %T in Global section)
  local vecscatter=$(grep "^VecScatterBegin " "$log" | awk '{print $11}')

  echo "${q_time:-n/a},${p_time:-n/a},${matmult:-n/a},${vecscatter:-n/a},${dofs:-n/a}"
}

echo "================================"
echo "Test 1: Problem Size Scaling"
echo "================================"
echo ""
echo "Config,DOFs,CPU_Quandary(s),CPU_PETSc(s),CPU_MatMult%,CPU_Scatter%,GPU_Quandary(s),GPU_PETSc(s),GPU_MatMult%,GPU_Scatter%,Speedup"

for nlevels in 4 16 32; do
  cpu_dir="${BENCH_DIR}/size_${nlevels}_cpu"
  gpu_dir="${BENCH_DIR}/size_${nlevels}_gpu"

  if [ ! -d "$cpu_dir" ] || [ ! -d "$gpu_dir" ]; then
    continue
  fi

  cpu_data=$(extract_from_run "$cpu_dir" "cpu")
  gpu_data=$(extract_from_run "$gpu_dir" "gpu")

  cpu_qt=$(echo "$cpu_data" | cut -d, -f1)
  cpu_pt=$(echo "$cpu_data" | cut -d, -f2)
  cpu_mm=$(echo "$cpu_data" | cut -d, -f3)
  cpu_vs=$(echo "$cpu_data" | cut -d, -f4)
  dofs=$(echo "$cpu_data" | cut -d, -f5)

  gpu_qt=$(echo "$gpu_data" | cut -d, -f1)
  gpu_pt=$(echo "$gpu_data" | cut -d, -f2)
  gpu_mm=$(echo "$gpu_data" | cut -d, -f3)
  gpu_vs=$(echo "$gpu_data" | cut -d, -f4)

  if [ "$cpu_qt" != "n/a" ] && [ "$gpu_qt" != "n/a" ]; then
    speedup=$(echo "scale=2; $cpu_qt / $gpu_qt" | bc)
  else
    speedup="n/a"
  fi

  echo "${nlevels}^4,${dofs},${cpu_qt},${cpu_pt},${cpu_mm},${cpu_vs},${gpu_qt},${gpu_pt},${gpu_mm},${gpu_vs},${speedup}x"
done

echo ""
echo "================================"
echo "Test 2: Strong Scaling (32^4)"
echo "================================"
echo ""
echo "Ranks,CPU_Quandary(s),CPU_PETSc(s),CPU_MatMult%,CPU_Scatter%,GPU_Quandary(s),GPU_PETSc(s),GPU_MatMult%,GPU_Scatter%"

for nranks in 1 4 8; do
  cpu_dir="${BENCH_DIR}/scaling_${nranks}ranks_cpu"
  gpu_dir="${BENCH_DIR}/scaling_${nranks}ranks_gpu"

  if [ ! -d "$cpu_dir" ]; then
    continue
  fi

  cpu_data=$(extract_from_run "$cpu_dir" "cpu")
  cpu_qt=$(echo "$cpu_data" | cut -d, -f1)
  cpu_pt=$(echo "$cpu_data" | cut -d, -f2)
  cpu_mm=$(echo "$cpu_data" | cut -d, -f3)
  cpu_vs=$(echo "$cpu_data" | cut -d, -f4)

  if [ -d "$gpu_dir" ]; then
    gpu_data=$(extract_from_run "$gpu_dir" "gpu")
    gpu_qt=$(echo "$gpu_data" | cut -d, -f1)
    gpu_pt=$(echo "$gpu_data" | cut -d, -f2)
    gpu_mm=$(echo "$gpu_data" | cut -d, -f3)
    gpu_vs=$(echo "$gpu_data" | cut -d, -f4)
  else
    gpu_qt="n/a"
    gpu_pt="n/a"
    gpu_mm="n/a"
    gpu_vs="n/a"
  fi

  echo "${nranks},${cpu_qt},${cpu_pt},${cpu_mm},${cpu_vs},${gpu_qt},${gpu_pt},${gpu_mm},${gpu_vs}"
done

echo ""
echo "================================"
echo "Key Metrics Guide"
echo "================================"
echo "Quandary time: Solve loop only (excludes setup)"
echo "PETSc time: Total including initialization"
echo "MatMult%: Compute time (higher = compute-bound)"
echo "Scatter%: Communication time (higher = communication-bound)"
echo ""
echo "Good GPU performance: MatMult% < 60%, Scatter% < 25%, Speedup > 5x"
