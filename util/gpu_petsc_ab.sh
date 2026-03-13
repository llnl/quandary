#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Build a small PETSc GPU/CPU matrix (Spack envs) and run Quandary with PETSc logging.

Defaults are set up for LLNL Tioga using radiuss-spack-configs.

Usage:
  bash util/gpu_petsc_ab.sh [options]

Options:
  --machine {tioga|tuolumne}     Select radiuss machine config (default: tioga)
  --repo PATH                   Path to Quandary repo root (default: cwd)
  --develop                     Use Spack develop mode for Quandary (default: on)
  --no-develop                  Disable Spack develop mode for Quandary
  --develop-path PATH           Path to local Quandary source for develop mode (default: --repo)
  --env-root PATH               Where to create envs/ (default: <repo>/envs-gpu-ab)
  --variants LIST               Which variants to build/run (default: cpu,kokkos,rocm)
                                LIST is comma-separated, e.g. "cpu,kokkos" or "kokkos"
  --run-only                    Do not modify/concretize/install envs; just run existing binaries
  --no-quiet                    Do not pass --quiet to Quandary (default: pass --quiet)
  --launcher "CMD"              Launcher prefix (default: flux run; adds --gpu-bind=closest if supported)
                                The script appends "-n <nprocs>" automatically.
  --no-mpich-gpu-support        Do not set MPICH_GPU_SUPPORT_ENABLED=1 for GPU variants
  --petsc-gpu-aware-mpi MODE    Control PETSc GPU-aware MPI behavior for GPU variants
                                MODE is {auto|on|off} (default: off for kokkos; auto otherwise)
  --petsc-debug                 Build PETSc with +debug (default: ~debug)
  --debugger "CMD"              Run Quandary under a debugger wrapper (disables /usr/bin/time)
                                Example: --debugger "gdb -batch -ex run -ex bt --args"
  --kokkos-vec-type TYPE        PETSc vec type for kokkos variant (default: kokkos)
  --kokkos-mat-type TYPE        PETSc mat type for kokkos variant (default: aijkokkos)
  --run-output-root PATH        Root directory for Quandary outputs (default: <env-root>/runs)
  --no-unique-output            Do not rewrite config output directory (default: unique per run under --run-output-root)
  --nprocs N                    MPI ranks (default: 8 on tioga, 4 on tuolumne)
  --cfg PATH                    Config file (default: tests/performance/configs/nlevels_32_32_32_32.toml)
  --llvm-amdgpu VER             llvm-amdgpu compiler version (default: 6.4.3)
  --hip VER                     hip/ROCm version (default: 6.4.3; used for GPU variants)
  --pin-hipblas                 Require hipblas@--hip (errors if unavailable; default: auto-detect if available)
  --no-pin-hipblas              Never pin hipblas version (default: auto-detect)
  --petsc-cxxstd N              PETSc C++ standard (adds PETSc cxxflags="-std=gnu++N"); default: unset
  --kokkos-cxxstd N             Kokkos C++ standard for PETSc+Kokkos builds (default: 17)
  --amdgpu-target GFX           AMDGPU target (default: gfx90a on tioga, gfx942 on tuolumne)
  --petsc SPEC                  PETSc version/constraint appended after "^petsc"
                                (default: "@3.24.4", e.g. "@3.24.4" or "@3.24.4:")
  --petsc-min VARS              Extra PETSc variant toggles (default: "~mmg~parmmg~saws~examples~ml~exodusii~zoltan")
  --with-test-deps              Build Quandary with "+test" (adds python/pip run deps; default: off)
  --no-install                  Only concretize; skip spack install
  --keep-logs                   Don’t overwrite logs; append timestamp
  -h, --help                    Show help

Outputs:
  <env-root>/cpu.log
  <env-root>/kokkos.log
  <env-root>/rocm.log
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

REPO="${PWD}"
MACHINE="tioga"
ENV_ROOT=""
DEVELOP="1"
DEVELOP_PATH=""
VARIANTS="cpu,kokkos,rocm"
LAUNCHER=""
LAUNCHER_SET="0"
NPROCS=""
CFG="tests/performance/configs/nlevels_32_32_32_32.toml"
LLVM_AMDGPU_VER="6.4.3"
HIP_VER="6.4.3"
AMDGPU_TARGET=""
PETSC_SPEC="@3.24.4"
PETSC_MIN="~mmg~parmmg~saws~examples~ml~exodusii~zoltan"
PETSC_DEBUG="~debug"
PETSC_CXXSTD=""
KOKKOS_CXXSTD="17"
QUANDARY_TEST_VARIANT="~test"
DO_INSTALL="1"
KEEP_LOGS="0"
RUN_ONLY="0"
QUANDARY_QUIET="1"
HIPBLAS_PIN_MODE="auto" # auto|pin|off
MPICH_GPU_SUPPORT="1"
PETSC_GPU_AWARE_MPI="auto" # auto|on|off
PETSC_GPU_AWARE_MPI_SET="0"
KOKKOS_VEC_TYPE="kokkos"
KOKKOS_MAT_TYPE="aijkokkos"
DEBUGGER=""
RUN_OUTPUT_ROOT=""
UNIQUE_OUTPUT="1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --machine) MACHINE="$2"; shift 2;;
    --repo) REPO="$2"; shift 2;;
    --develop) DEVELOP="1"; shift 1;;
    --no-develop) DEVELOP="0"; shift 1;;
    --develop-path) DEVELOP_PATH="$2"; shift 2;;
    --env-root) ENV_ROOT="$2"; shift 2;;
    --variants) VARIANTS="$2"; shift 2;;
    --run-only) RUN_ONLY="1"; shift 1;;
    --no-quiet) QUANDARY_QUIET="0"; shift 1;;
    --launcher) LAUNCHER="$2"; LAUNCHER_SET="1"; shift 2;;
    --no-mpich-gpu-support) MPICH_GPU_SUPPORT="0"; shift 1;;
    --petsc-gpu-aware-mpi) PETSC_GPU_AWARE_MPI="$2"; PETSC_GPU_AWARE_MPI_SET="1"; shift 2;;
    --petsc-debug) PETSC_DEBUG="+debug"; shift 1;;
    --debugger) DEBUGGER="$2"; shift 2;;
    --kokkos-vec-type) KOKKOS_VEC_TYPE="$2"; shift 2;;
    --kokkos-mat-type) KOKKOS_MAT_TYPE="$2"; shift 2;;
    --run-output-root) RUN_OUTPUT_ROOT="$2"; shift 2;;
    --no-unique-output) UNIQUE_OUTPUT="0"; shift 1;;
    --nprocs) NPROCS="$2"; shift 2;;
    --cfg) CFG="$2"; shift 2;;
    --llvm-amdgpu) LLVM_AMDGPU_VER="$2"; shift 2;;
    --hip) HIP_VER="$2"; shift 2;;
    --pin-hipblas) HIPBLAS_PIN_MODE="pin"; shift 1;;
    --no-pin-hipblas) HIPBLAS_PIN_MODE="off"; shift 1;;
    --petsc-cxxstd) PETSC_CXXSTD="$2"; shift 2;;
    --kokkos-cxxstd) KOKKOS_CXXSTD="$2"; shift 2;;
    --amdgpu-target) AMDGPU_TARGET="$2"; shift 2;;
    --petsc) PETSC_SPEC="$2"; shift 2;;
    --petsc-min) PETSC_MIN="$2"; shift 2;;
    --with-test-deps) QUANDARY_TEST_VARIANT="+test"; shift 1;;
    --no-install) DO_INSTALL="0"; shift 1;;
    --keep-logs) KEEP_LOGS="1"; shift 1;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

command -v spack >/dev/null 2>&1 || die "'spack' not found in PATH"

case "$PETSC_GPU_AWARE_MPI" in
  auto|on|off) ;;
  *) die "--petsc-gpu-aware-mpi must be one of: auto, on, off";;
esac

cd "$REPO"
[[ -f "CMakeLists.txt" ]] || die "--repo must point to the Quandary repo root (missing CMakeLists.txt)"
: "${DEVELOP_PATH:=$REPO}"

default_launcher() {
  local launcher="flux run --gpus-per-task=1"
  if command -v flux >/dev/null 2>&1; then
    if flux run --help 2>&1 | grep -q -- '--gpu-bind'; then
      launcher+=" --gpu-bind=closest"
    fi
  fi
  echo "$launcher"
}

if [[ -z "$LAUNCHER" ]]; then
  LAUNCHER="$(default_launcher)"
fi

RAD="$REPO/.ci-scripts/radiuss-spack-configs/toss_4_x86_64_ib_cray"
[[ -d "$RAD" ]] || die "radiuss-spack-configs not found at: $RAD"

case "$MACHINE" in
  tioga)
    : "${NPROCS:=8}"
    : "${AMDGPU_TARGET:=gfx90a}"
    MACHINE_PACKAGES="$RAD/tioga/packages.yaml"
    ;;
  tuolumne)
    : "${NPROCS:=4}"
    : "${AMDGPU_TARGET:=gfx942}"
    MACHINE_PACKAGES="$RAD/tuolumne/packages.yaml"
    ;;
  *)
    die "--machine must be tioga or tuolumne"
    ;;
esac

[[ -f "$MACHINE_PACKAGES" ]] || die "Missing machine packages.yaml: $MACHINE_PACKAGES"
[[ -f "$CFG" ]] || die "Config not found: $CFG"

if [[ -z "$ENV_ROOT" ]]; then
  ENV_ROOT="$REPO/envs-gpu-ab"
fi
if [[ -z "$RUN_OUTPUT_ROOT" ]]; then
  RUN_OUTPUT_ROOT="$ENV_ROOT/runs"
fi

if [[ "$RUN_ONLY" == "1" ]]; then
  DO_INSTALL="0"
fi

normalize_variants() {
  local raw="$1"
  raw="${raw// /}"
  echo "$raw"
}

variant_selected() {
  local needle="$1"
  local list="$2"
  local IFS=,
  local v
  for v in $list; do
    [[ "$v" == "$needle" ]] && return 0
  done
  return 1
}

validate_variants() {
  local list="$1"
  local IFS=,
  local v
  for v in $list; do
    case "$v" in
      cpu|kokkos|rocm) ;;
      *) die "Invalid --variants entry: '$v' (allowed: cpu,kokkos,rocm)";;
    esac
  done
}

VARIANTS="$(normalize_variants "$VARIANTS")"
[[ -n "$VARIANTS" ]] || die "--variants cannot be empty"
validate_variants "$VARIANTS"

mkdir -p "$ENV_ROOT"

if [[ "$PETSC_GPU_AWARE_MPI_SET" == "0" ]] && variant_selected kokkos "$VARIANTS"; then
  # The PETSc+Kokkos path often triggers MPI communication with device-resident
  # buffers; default to PETSc staging (off) unless the user explicitly opts in.
  PETSC_GPU_AWARE_MPI="off"
fi

write_env_yaml() {
  local env_dir="$1"
  cat > "${env_dir}/spack.yaml" <<EOF
spack:
  include:
    - ${RAD}/config.yaml
    - ${MACHINE_PACKAGES}
  concretizer:
    unify: true
  view: true
EOF
  if [[ "$DEVELOP" == "1" ]]; then
    cat >> "${env_dir}/spack.yaml" <<EOF
  develop:
    quandary:
      path: ${DEVELOP_PATH}
      spec: quandary@main
EOF
  fi
  cat >> "${env_dir}/spack.yaml" <<EOF
  specs: []
EOF
}

if [[ "$RUN_ONLY" != "1" ]]; then
  if variant_selected cpu "$VARIANTS"; then
    mkdir -p "$ENV_ROOT/cpu"
    write_env_yaml "$ENV_ROOT/cpu"
  fi
  if variant_selected kokkos "$VARIANTS"; then
    mkdir -p "$ENV_ROOT/kokkos"
    write_env_yaml "$ENV_ROOT/kokkos"
  fi
  if variant_selected rocm "$VARIANTS"; then
    mkdir -p "$ENV_ROOT/rocm"
    write_env_yaml "$ENV_ROOT/rocm"
  fi
else
  if variant_selected cpu "$VARIANTS"; then [[ -d "$ENV_ROOT/cpu" ]] || die "Missing env dir for --run-only: $ENV_ROOT/cpu"; fi
  if variant_selected kokkos "$VARIANTS"; then [[ -d "$ENV_ROOT/kokkos" ]] || die "Missing env dir for --run-only: $ENV_ROOT/kokkos"; fi
  if variant_selected rocm "$VARIANTS"; then [[ -d "$ENV_ROOT/rocm" ]] || die "Missing env dir for --run-only: $ENV_ROOT/rocm"; fi
fi

echo "=== Env status ==="
if variant_selected cpu "$VARIANTS"; then spack -e "$ENV_ROOT/cpu" env status || true; fi
if variant_selected kokkos "$VARIANTS"; then spack -e "$ENV_ROOT/kokkos" env status || true; fi
if variant_selected rocm "$VARIANTS"; then spack -e "$ENV_ROOT/rocm" env status || true; fi

COMPILER_SPEC="%llvm-amdgpu@=${LLVM_AMDGPU_VER}"

hipblas_version_available() {
  local ver="$1"
  # If Spack doesn't know about hipblas, we'll find out later during concretize.
  # Here we only guard against pinning an unavailable version.
  spack versions -s hipblas 2>/dev/null | grep -q -E "(^|[[:space:]])${ver}([[:space:]]|$)"
}

HIPBLAS_SPEC="^hipblas"
case "$HIPBLAS_PIN_MODE" in
  off)
    ;;
  pin)
    if hipblas_version_available "${HIP_VER}"; then
      HIPBLAS_SPEC="^hipblas@${HIP_VER}"
    else
      die "hipblas@${HIP_VER} is not available in this Spack; omit --pin-hipblas (auto mode) or use --no-pin-hipblas"
    fi
    ;;
  auto)
    if hipblas_version_available "${HIP_VER}"; then
      HIPBLAS_SPEC="^hipblas@${HIP_VER}"
    fi
    ;;
  *)
    die "Internal error: invalid HIPBLAS_PIN_MODE='${HIPBLAS_PIN_MODE}'"
    ;;
esac

# hipblas-common is a separate package on some stacks; keeping it aligned avoids
# mixing ROCm component versions (e.g. hipblas@6.4.x with hipblas-common@7.x).
# Do not apply a compiler constraint here; on systems where hipblas-common is
# provided externally by ROCm, a compiler constraint can prevent matching it.
HIPBLAS_COMMON_SPEC="^hipblas-common@${HIP_VER}"

# For GPU variants, keep ROCm deps consistent. Do not force a compiler on ROCm
# externals like hipblas-common; some site stacks provide them prebuilt.
ROCM_DEPS="^hip@${HIP_VER} ${HIPBLAS_SPEC} ${HIPBLAS_COMMON_SPEC}"

# Make PETSc explicitly use the same compiler as Quandary. Place the compiler
# constraint at the end so PETSc variants don't get parsed as compiler variants.
PETSC_CPU="^petsc${PETSC_SPEC}~rocm~kokkos ${PETSC_MIN} ${COMPILER_SPEC}"
PETSC_KOKKOS="^petsc${PETSC_SPEC}+rocm+kokkos amdgpu_target=${AMDGPU_TARGET} ${PETSC_MIN} ${COMPILER_SPEC}"
PETSC_ROCM="^petsc${PETSC_SPEC}+rocm~kokkos amdgpu_target=${AMDGPU_TARGET} ${PETSC_MIN} ${COMPILER_SPEC}"

PETSC_CPU="${PETSC_CPU/^petsc${PETSC_SPEC}/^petsc${PETSC_SPEC}${PETSC_DEBUG}}"
PETSC_KOKKOS="${PETSC_KOKKOS/^petsc${PETSC_SPEC}/^petsc${PETSC_SPEC}${PETSC_DEBUG}}"
PETSC_ROCM="${PETSC_ROCM/^petsc${PETSC_SPEC}/^petsc${PETSC_SPEC}${PETSC_DEBUG}}"

if [[ -n "$PETSC_CXXSTD" ]]; then
  # Spack's PETSc package doesn't necessarily expose a cxxstd variant; use flags.
  # Keep this before the compiler constraint so it applies to PETSc, not the compiler spec.
  PETSC_KOKKOS="${PETSC_KOKKOS%${COMPILER_SPEC}} cxxflags=-std=gnu++${PETSC_CXXSTD} ${COMPILER_SPEC}"
fi

KOKKOS_CONSTRAINT="^kokkos cxxstd=${KOKKOS_CXXSTD}"

if [[ "$RUN_ONLY" != "1" ]]; then
  echo "=== (Re)setting specs ==="
  if variant_selected cpu "$VARIANTS"; then
    spack -e "$ENV_ROOT/cpu" rm -y quandary >/dev/null 2>&1 || true
    # CPU baseline: do not mention hip or amdgpu_target; keep rocm disabled at the Quandary level.
    spack -e "$ENV_ROOT/cpu" add "quandary@main${QUANDARY_TEST_VARIANT}~rocm ${COMPILER_SPEC} ^cray-mpich${COMPILER_SPEC} ${PETSC_CPU}"
  fi
  if variant_selected kokkos "$VARIANTS"; then
    spack -e "$ENV_ROOT/kokkos" rm -y quandary >/dev/null 2>&1 || true
    # GPU variant: enable ROCm at the Quandary level to keep amdgpu_target consistent for PETSc deps.
    # Kokkos is selected on PETSc (Quandary has no +kokkos variant).
    # Keep Kokkos itself on C++17; PETSc's Kokkos package currently rejects C++20.
    spack -e "$ENV_ROOT/kokkos" add "quandary@main${QUANDARY_TEST_VARIANT}+rocm amdgpu_target=${AMDGPU_TARGET} ${COMPILER_SPEC} ^cray-mpich${COMPILER_SPEC} ${ROCM_DEPS} ${PETSC_KOKKOS} ${KOKKOS_CONSTRAINT}"
  fi
  if variant_selected rocm "$VARIANTS"; then
    spack -e "$ENV_ROOT/rocm" rm -y quandary >/dev/null 2>&1 || true
    # GPU variant: enable ROCm at the Quandary level to keep amdgpu_target consistent for PETSc deps.
    spack -e "$ENV_ROOT/rocm" add "quandary@main${QUANDARY_TEST_VARIANT}+rocm amdgpu_target=${AMDGPU_TARGET} ${COMPILER_SPEC} ^cray-mpich${COMPILER_SPEC} ${ROCM_DEPS} ${PETSC_ROCM}"
  fi

  echo "=== Concretize ==="
  if variant_selected cpu "$VARIANTS"; then spack -e "$ENV_ROOT/cpu" concretize -f; fi
  if variant_selected kokkos "$VARIANTS"; then spack -e "$ENV_ROOT/kokkos" concretize -f; fi
  if variant_selected rocm "$VARIANTS"; then spack -e "$ENV_ROOT/rocm" concretize -f; fi
else
  echo "=== Run only (--run-only): skipping spec reset/concretize/install ==="
fi

if [[ "$DO_INSTALL" == "1" ]]; then
  echo "=== Install ==="
  if variant_selected cpu "$VARIANTS"; then spack -e "$ENV_ROOT/cpu" install; fi
  if variant_selected kokkos "$VARIANTS"; then spack -e "$ENV_ROOT/kokkos" install; fi
  if variant_selected rocm "$VARIANTS"; then spack -e "$ENV_ROOT/rocm" install; fi
else
  echo "=== Skipping install (--no-install) ==="
fi

declare -A BIN_PATHS
declare -A PETSC_OPTS
declare -A LOG_PATHS

PETSC_OPTS[cpu]='-log_view -log_summary'
PETSC_OPTS[kokkos]="-vec_type ${KOKKOS_VEC_TYPE} -mat_type ${KOKKOS_MAT_TYPE} -log_view -log_summary"
PETSC_OPTS[rocm]='-vec_type hip -mat_type aijhipsparse -log_view -log_summary'

case "$PETSC_GPU_AWARE_MPI" in
  on)
    PETSC_OPTS[kokkos]="${PETSC_OPTS[kokkos]} -use_gpu_aware_mpi 1"
    PETSC_OPTS[rocm]="${PETSC_OPTS[rocm]} -use_gpu_aware_mpi 1"
    ;;
  off)
    PETSC_OPTS[kokkos]="${PETSC_OPTS[kokkos]} -use_gpu_aware_mpi 0"
    PETSC_OPTS[rocm]="${PETSC_OPTS[rocm]} -use_gpu_aware_mpi 0"
    ;;
  auto) ;;
esac

if variant_selected cpu "$VARIANTS"; then
  BIN_PATHS[cpu]="$(spack -e "$ENV_ROOT/cpu" location -i quandary)/bin/quandary"
  [[ -x "${BIN_PATHS[cpu]}" ]] || die "Binary not found/executable: ${BIN_PATHS[cpu]}"
fi
if variant_selected kokkos "$VARIANTS"; then
  BIN_PATHS[kokkos]="$(spack -e "$ENV_ROOT/kokkos" location -i quandary)/bin/quandary"
  [[ -x "${BIN_PATHS[kokkos]}" ]] || die "Binary not found/executable: ${BIN_PATHS[kokkos]}"
fi
if variant_selected rocm "$VARIANTS"; then
  BIN_PATHS[rocm]="$(spack -e "$ENV_ROOT/rocm" location -i quandary)/bin/quandary"
  [[ -x "${BIN_PATHS[rocm]}" ]] || die "Binary not found/executable: ${BIN_PATHS[rocm]}"
fi

logname() {
  local base="$1"
  if [[ "$KEEP_LOGS" == "1" ]]; then
    echo "${base%.*}_$(date +%Y%m%d_%H%M%S).log"
  else
    echo "$base"
  fi
}

if variant_selected cpu "$VARIANTS"; then LOG_PATHS[cpu]="${ENV_ROOT}/$(logname cpu.log)"; fi
if variant_selected kokkos "$VARIANTS"; then LOG_PATHS[kokkos]="${ENV_ROOT}/$(logname kokkos.log)"; fi
if variant_selected rocm "$VARIANTS"; then LOG_PATHS[rocm]="${ENV_ROOT}/$(logname rocm.log)"; fi

MPI_PREFIX="${LAUNCHER} -n ${NPROCS}"

echo "=== Run settings ==="
echo "machine=$MACHINE nprocs=$NPROCS"
echo "launcher='$LAUNCHER'"
echo "variants=$VARIANTS"
echo "cfg=$CFG"
echo "env_root=$ENV_ROOT"
if [[ "$UNIQUE_OUTPUT" == "1" ]]; then echo "run_output_root=$RUN_OUTPUT_ROOT"; fi
if [[ -n "$PETSC_CXXSTD" ]]; then echo "petsc_cxxstd=$PETSC_CXXSTD"; fi
if variant_selected kokkos "$VARIANTS"; then echo "kokkos_cxxstd=$KOKKOS_CXXSTD"; fi
if [[ "$PETSC_GPU_AWARE_MPI" != "auto" ]]; then echo "petsc_gpu_aware_mpi=$PETSC_GPU_AWARE_MPI"; fi
echo

set -x

run_variant() {
  local v="$1"
  local cfg_path="$CFG"
  if [[ "$UNIQUE_OUTPUT" == "1" ]]; then
    local ts
    ts="$(date +%Y%m%d_%H%M%S)"
    local tag="${v}"
    if [[ "$v" == "kokkos" ]]; then
      tag+="_vec-${KOKKOS_VEC_TYPE}_mat-${KOKKOS_MAT_TYPE}"
    elif [[ "$v" == "rocm" ]]; then
      tag+="_vec-hip_mat-aijhipsparse"
    fi
    tag="$(echo "$tag" | tr -c 'A-Za-z0-9._-' '_' )"
    local outdir="${RUN_OUTPUT_ROOT}/${tag}_${ts}"
    mkdir -p "$outdir"
    cfg_path="${outdir}/config.toml"
    awk -v out="$outdir" '
      BEGIN { in_output = 0 }
      /^\[/ { in_output = ($0 == "[output]") }
      in_output && $0 ~ /^[[:space:]]*directory[[:space:]]*=/ { print "directory = \"" out "\""; next }
      { print }
    ' "$CFG" > "$cfg_path"
  fi
  local -a cmd_prefix
  if [[ "$MPICH_GPU_SUPPORT" == "1" ]] && [[ "$v" != "cpu" ]]; then
    # Cray MPICH typically requires this for GPU-aware MPI. Without it, PETSc may
    # pass device pointers through MPI and trigger invalid memory accesses.
    cmd_prefix=(env MPICH_GPU_SUPPORT_ENABLED=1)
  fi
  if [[ -n "$DEBUGGER" ]]; then
    # Note: for MPI runs, this will launch a debugger per rank. Prefer --nprocs 1.
    "${cmd_prefix[@]}" ${MPI_PREFIX} ${DEBUGGER} "${BIN_PATHS[$v]}" "$cfg_path" \
      --petsc-options "${PETSC_OPTS[$v]}" 2>&1 | tee "${LOG_PATHS[$v]}"
    return 0
  fi
  if [[ "$QUANDARY_QUIET" == "1" ]]; then
    "${cmd_prefix[@]}" ${MPI_PREFIX} /usr/bin/time -p "${BIN_PATHS[$v]}" "$cfg_path" --quiet \
      --petsc-options "${PETSC_OPTS[$v]}" 2>&1 | tee "${LOG_PATHS[$v]}"
  else
    "${cmd_prefix[@]}" ${MPI_PREFIX} /usr/bin/time -p "${BIN_PATHS[$v]}" "$cfg_path" \
      --petsc-options "${PETSC_OPTS[$v]}" 2>&1 | tee "${LOG_PATHS[$v]}"
  fi
}

if variant_selected cpu "$VARIANTS"; then run_variant cpu; fi
if variant_selected kokkos "$VARIANTS"; then run_variant kokkos; fi
if variant_selected rocm "$VARIANTS"; then run_variant rocm; fi

set +x

echo "Done. Logs:"
if variant_selected cpu "$VARIANTS"; then
  echo "  ${LOG_PATHS[cpu]}"
fi
if variant_selected kokkos "$VARIANTS"; then
  echo "  ${LOG_PATHS[kokkos]}"
fi
if variant_selected rocm "$VARIANTS"; then
  echo "  ${LOG_PATHS[rocm]}"
fi
