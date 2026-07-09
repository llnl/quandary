"""Shared pytest fixtures for Python tests."""

import os
from pathlib import Path
import shlex
import shutil
from types import SimpleNamespace

import pytest


def _normalize_launcher_value(value: str | None):
    """Convert CMake-style command lists into a shell-style launcher string."""
    if not value:
        return None
    return value.replace(";", " ").strip() or None


@pytest.fixture(scope="session", autouse=True)
def quandary_on_path():
    """Make the repo-root quandary binary discoverable without manual PATH setup."""
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    os.environ["PATH"] = repo_root + os.pathsep + os.environ.get("PATH", "")


def _parse_cmake_cache(cache_path: Path):
    """Read MPI launcher settings from a configured CMake cache if available."""
    if not cache_path.is_file():
        return None, None

    mpi_exec = None
    nproc_flag = None
    with cache_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("MPIEXEC_EXECUTABLE:") or line.startswith("MPIEXEC:PATH=") or line.startswith("MPIEXEC:FILEPATH="):
                mpi_exec = _normalize_launcher_value(line.split("=", 1)[1].strip())
            elif line.startswith("MPIEXEC_NUMPROC_FLAG:"):
                nproc_flag = line.split("=", 1)[1].strip()
    return mpi_exec or None, nproc_flag or None


def _parse_host_config(host_config_path: Path):
    """Read MPI launcher settings from a BLT host-config if available."""
    if not host_config_path.is_file():
        return None

    for raw_line in host_config_path.read_text().splitlines():
        line = raw_line.strip()
        if not line.startswith("set("):
            continue
        if "MPIEXEC_EXECUTABLE" in line or line.startswith("set(MPIEXEC "):
            parts = line.split('"')
            if len(parts) >= 2 and parts[1].strip():
                return _normalize_launcher_value(parts[1].strip())
    return None


def _infer_nproc_flag(mpi_exec: str) -> str:
    """Infer the process-count flag from the launcher executable name."""
    launcher_name = Path(shlex.split(_normalize_launcher_value(mpi_exec))[0]).name
    if launcher_name in {"mpirun", "mpiexec", "orterun"}:
        return "-np"
    return "-n"


def _launcher_name(mpi_exec: str) -> str:
    """Return the basename of the launcher executable."""
    return Path(shlex.split(_normalize_launcher_value(mpi_exec))[0]).name


def _has_flux_allocation() -> bool:
    """Detect whether the current process is already running inside a Flux allocation."""
    return any(
        os.environ.get(name)
        for name in ("FLUX_URI", "FLUX_JOB_ID", "FLUX_HANDLE")
    )


def _resolve_host_config_path(repo_root: Path):
    """Locate a generated or explicitly selected host-config file."""
    host_config = os.environ.get("HOST_CONFIG", "").strip()
    if host_config:
        host_config_path = Path(host_config)
        if not host_config_path.is_absolute():
            host_config_path = repo_root / host_config
        if host_config_path.is_file():
            return host_config_path

    top_level_configs = [
        path for path in repo_root.glob("*.cmake")
        if path.name != "CMakeLists.txt"
    ]
    if len(top_level_configs) == 1:
        return top_level_configs[0]

    return None


def resolve_mpi_launcher(request):
    """Resolve MPI launcher settings from pytest options, build cache, or host-config."""
    requested_exec = _normalize_launcher_value(request.config.getoption("--mpi-exec").strip()) or "mpirun"
    mpi_opt = request.config.getoption("--mpi-opt").strip()

    repo_root = Path(__file__).resolve().parents[2]
    mpi_exec = requested_exec
    nproc_flag = None

    if requested_exec == "mpirun":
        cache_exec, cache_flag = _parse_cmake_cache(repo_root / "build" / "CMakeCache.txt")
        if cache_exec:
            mpi_exec = cache_exec
            nproc_flag = cache_flag
        else:
            host_config_path = _resolve_host_config_path(repo_root)
            if host_config_path is not None:
                host_exec = _parse_host_config(host_config_path)
                if host_exec:
                    mpi_exec = host_exec

        if _launcher_name(mpi_exec) == "flux" and not _has_flux_allocation():
            srun_path = shutil.which("srun")
            if srun_path and (os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CLUSTER_NAME")):
                mpi_exec = srun_path
                nproc_flag = None

        if mpi_exec == "mpirun" and (os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CLUSTER_NAME")):
            srun_path = shutil.which("srun")
            if srun_path:
                mpi_exec = srun_path

    if nproc_flag is None:
        nproc_flag = _infer_nproc_flag(mpi_exec)

    mpi_exec = _normalize_launcher_value(mpi_exec) or "mpirun"
    mpi_launcher = " ".join(part for part in [mpi_exec, mpi_opt] if part)
    return SimpleNamespace(exec=mpi_launcher, nproc_flag=nproc_flag)


@pytest.fixture
def mpi_exec(request):
    """Get MPI executor and options from pytest options."""
    launcher = resolve_mpi_launcher(request)
    return f"{launcher.exec} {launcher.nproc_flag}"


@pytest.fixture
def mpi_launcher(request):
    """Resolved MPI launcher for tests that pass exec and process-count flag separately."""
    return resolve_mpi_launcher(request)
