#!/usr/bin/env bash
# Extract key metrics from benchmark logs
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: extract_results.sh [--table|--csv] <benchmark_dir>

Default output is a human-aligned table (uses `column -t -s,` when available).

Options:
  --table   Aligned table output (default)
  --csv     Raw CSV output (machine-friendly)
  -h,--help Show this help
EOF
}

FORMAT="table"
BENCH_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --table) FORMAT="table"; shift;;
    --csv) FORMAT="csv"; shift;;
    -h|--help) usage; exit 0;;
    *)
      if [[ -z "$BENCH_DIR" ]]; then
        BENCH_DIR="$1"
        shift
      else
        echo "ERROR: unexpected argument: $1" >&2
        usage >&2
        exit 1
      fi
      ;;
  esac
done

if [[ -z "$BENCH_DIR" ]]; then
  usage >&2
  exit 1
fi

maybe_column() {
  if [[ "$FORMAT" == "table" ]] && command -v column >/dev/null 2>&1; then
    column -t -s,
  else
    cat
  fi
}

is_number() {
  # Accepts integers, decimals, and scientific notation (e.g. 1.23e-04).
  [[ "${1:-}" =~ ^[+-]?[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$ ]] || [[ "${1:-}" =~ ^[+-]?[0-9]*[.][0-9]+([eE][+-]?[0-9]+)?$ ]]
}

calc_speedup() {
  local cpu="$1"
  local gpu="$2"

  if ! is_number "$cpu" || ! is_number "$gpu"; then
    echo "n/a"
    return
  fi

  awk -v cpu="$cpu" -v gpu="$gpu" 'BEGIN { if (gpu == 0) { print "n/a" } else { printf "%.2f", cpu / gpu } }'
}

calc_from_pct() {
  local total="$1"
  local pct="$2"

  if ! is_number "$total" || ! is_number "$pct"; then
    echo "n/a"
    return
  fi

  awk -v total="$total" -v pct="$pct" 'BEGIN { printf "%.6g", total * (pct / 100.0) }'
}

extract_used_time() {
  local log="$1"
  awk '/Used Time:/ {print $3; exit}' "$log" 2>/dev/null || true
}

extract_petsc_total_time() {
  local log="$1"
  grep "^Time (sec):" "$log" 2>/dev/null | awk '{print $3}' | head -n 1
}

extract_event_time() {
  local log="$1"
  local event="$2"
  grep "^${event} " "$log" 2>/dev/null | awk '{print $4}' | head -n 1
}

extract_event_pct() {
  local log="$1"
  local event="$2"
  grep "^${event} " "$log" 2>/dev/null | awk '{print $11}' | head -n 1
}

extract_event_time_or_estimate() {
  local log="$1"
  local event="$2"
  local petsc_total_time="$3"

  local t
  t=$(extract_event_time "$log" "$event")
  if is_number "$t"; then
    echo "$t"
    return
  fi

  local pct
  pct=$(extract_event_pct "$log" "$event")
  if is_number "$pct"; then
    calc_from_pct "$petsc_total_time" "$pct"
    return
  fi

  echo "n/a"
}

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

  # Quandary time from timing.dat (if available); otherwise fall back to the "Used Time" line in the log.
  local q_time="n/a"
  if [ -n "$run_subdir" ] && [ -f "$run_subdir/timing.dat" ]; then
    q_time=$(awk '{print $2}' "$run_subdir/timing.dat")
  elif [ -f "$log" ]; then
    q_time=$(extract_used_time "$log")
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
  local p_time
  p_time=$(extract_petsc_total_time "$log")

  # MatMult %T (column 11 is %T in Global section)
  local matmult=$(grep "^MatMult " "$log" | awk '{print $11}')

  # VecScatterBegin %T (column 11 is %T in Global section)
  local vecscatter=$(grep "^VecScatterBegin " "$log" | awk '{print $11}')

  echo "${q_time:-n/a},${p_time:-n/a},${matmult:-n/a},${vecscatter:-n/a},${dofs:-n/a}"
}

# If the provided directory is itself a benchmark result directory containing cpu.log/gpu.log,
# emit a compact summary for just that run and exit.
if [[ -f "${BENCH_DIR}/cpu.log" || -f "${BENCH_DIR}/gpu.log" ]]; then
  cpu_log="${BENCH_DIR}/cpu.log"
  gpu_log="${BENCH_DIR}/gpu.log"

  cpu_overall="n/a"
  gpu_overall="n/a"
  cpu_petsc="n/a"
  gpu_petsc="n/a"

  if [[ -f "$cpu_log" ]]; then
    cpu_overall=$(extract_used_time "$cpu_log")
    cpu_petsc=$(extract_petsc_total_time "$cpu_log")
  fi
  if [[ -f "$gpu_log" ]]; then
    gpu_overall=$(extract_used_time "$gpu_log")
    gpu_petsc=$(extract_petsc_total_time "$gpu_log")
  fi

  cpu_matmult=$(extract_event_time_or_estimate "$cpu_log" "MatMult" "$cpu_petsc")
  gpu_matmult=$(extract_event_time_or_estimate "$gpu_log" "MatMult" "$gpu_petsc")

  cpu_vs_begin=$(extract_event_time_or_estimate "$cpu_log" "VecScatterBegin" "$cpu_petsc")
  cpu_vs_end=$(extract_event_time_or_estimate "$cpu_log" "VecScatterEnd" "$cpu_petsc")
  gpu_vs_begin=$(extract_event_time_or_estimate "$gpu_log" "VecScatterBegin" "$gpu_petsc")
  gpu_vs_end=$(extract_event_time_or_estimate "$gpu_log" "VecScatterEnd" "$gpu_petsc")

  cpu_vecscatter="n/a"
  gpu_vecscatter="n/a"
  if is_number "$cpu_vs_begin" && is_number "$cpu_vs_end"; then
    cpu_vecscatter=$(awk -v a="$cpu_vs_begin" -v b="$cpu_vs_end" 'BEGIN{printf "%.6g", a+b}')
  fi
  if is_number "$gpu_vs_begin" && is_number "$gpu_vs_end"; then
    gpu_vecscatter=$(awk -v a="$gpu_vs_begin" -v b="$gpu_vs_end" 'BEGIN{printf "%.6g", a+b}')
  fi

  echo "================================"
  echo "Single Run Summary"
  echo "================================"
  echo ""
  {
    echo "Metric,CPU(s),GPU(s),Speedup"
    s=$(calc_speedup "$cpu_overall" "$gpu_overall"); [[ "$s" != "n/a" ]] && s="${s}x"; echo "Overall,${cpu_overall:-n/a},${gpu_overall:-n/a},${s}"
    s=$(calc_speedup "$cpu_petsc" "$gpu_petsc"); [[ "$s" != "n/a" ]] && s="${s}x"; echo "PETSc Total,${cpu_petsc:-n/a},${gpu_petsc:-n/a},${s}"
    s=$(calc_speedup "$cpu_matmult" "$gpu_matmult"); [[ "$s" != "n/a" ]] && s="${s}x"; echo "MatMult,${cpu_matmult:-n/a},${gpu_matmult:-n/a},${s}"
    s=$(calc_speedup "$cpu_vecscatter" "$gpu_vecscatter"); [[ "$s" != "n/a" ]] && s="${s}x"; echo "VecScatter (Begin+End),${cpu_vecscatter:-n/a},${gpu_vecscatter:-n/a},${s}"
  } | maybe_column

  echo ""
  {
    echo "Detail,CPU(s),GPU(s)"
    echo "VecScatterBegin,${cpu_vs_begin:-n/a},${gpu_vs_begin:-n/a}"
    echo "VecScatterEnd,${cpu_vs_end:-n/a},${gpu_vs_end:-n/a}"
  } | maybe_column

  exit 0
fi

echo "================================"
echo "Test 1: Problem Size Scaling"
echo "================================"
echo ""
{
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

    speedup=$(calc_speedup "$cpu_qt" "$gpu_qt")
    if [ "$speedup" != "n/a" ]; then
      speedup="${speedup}x"
    fi

    echo "${nlevels}^4,${dofs},${cpu_qt},${cpu_pt},${cpu_mm},${cpu_vs},${gpu_qt},${gpu_pt},${gpu_mm},${gpu_vs},${speedup}"
  done
} | maybe_column

echo ""
echo "================================"
echo "Test 2: Strong Scaling (32^4)"
echo "================================"
echo ""
{
  echo "Ranks,CPU_Quandary(s),CPU_PETSc(s),CPU_MatMult%,CPU_Scatter%,GPU_Quandary(s),GPU_PETSc(s),GPU_MatMult%,GPU_Scatter%"

  # Find all scaling directories and extract rank numbers
  SCALING_DIRS=$(find "$BENCH_DIR" -maxdepth 1 -type d -name "scaling_*ranks_cpu" | sort -V)
  RANK_NUMS=()
  for dir in $SCALING_DIRS; do
    nranks=$(basename "$dir" | sed 's/scaling_\([0-9]*\)ranks_cpu/\1/')
    RANK_NUMS+=($nranks)
  done

  for nranks in "${RANK_NUMS[@]}"; do
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
} | maybe_column

echo ""
echo "================================"
echo "Test 3: Fixed GPUs (32^4)"
echo "================================"
echo ""
{
  echo "Ranks,GPUs,GPU_Quandary(s),GPU_PETSc(s),GPU_MatMult%,GPU_Scatter%"

  FIXED_DIRS=$(find "$BENCH_DIR" -maxdepth 1 -type d -name "fixed_*gpus_*ranks_gpu" | sort -V)
  for dir in $FIXED_DIRS; do
    base=$(basename "$dir")
    gpus=$(echo "$base" | sed -n 's/^fixed_\([0-9]\+\)gpus_.*$/\1/p')
    nranks=$(echo "$base" | sed -n 's/^fixed_[0-9]\+gpus_\([0-9]\+\)ranks_gpu$/\1/p')

    gpu_data=$(extract_from_run "$dir" "gpu")
    gpu_qt=$(echo "$gpu_data" | cut -d, -f1)
    gpu_pt=$(echo "$gpu_data" | cut -d, -f2)
    gpu_mm=$(echo "$gpu_data" | cut -d, -f3)
    gpu_vs=$(echo "$gpu_data" | cut -d, -f4)

    echo "${nranks:-n/a},${gpus:-n/a},${gpu_qt},${gpu_pt},${gpu_mm},${gpu_vs}"
  done
} | maybe_column

echo ""
echo "================================"
echo "Key Metrics Guide"
echo "================================"
echo "Quandary time: Solve loop only (excludes setup)"
echo "PETSc time: Total including initialization"
echo "MatMult%: Compute time (higher = compute-bound)"
echo "Scatter%: Communication time (higher = communication-bound)"
