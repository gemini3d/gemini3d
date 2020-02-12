import subprocess
import os
import shutil
from pathlib import Path
import math
import typing as T
try:
    import psutil
except ImportError:
    psutil = None
    # pip install psutil will improve CPU utilization.

from .raw import get_simsize as get_simsize_raw
from .hdf import get_simsize as get_simsize_h5

git = shutil.which("git")

Pathlike = T.Union[str, Path]

__all__ = ["gitrev", "get_cpu_count", "get_mpi_count", "get_simsize"]


def gitrev() -> str:
    if not git:
        return ""

    return subprocess.check_output([git, "rev-parse", "--short", "HEAD"], universal_newlines=True).strip()


def get_cpu_count(force: int = None) -> int:
    if force:
        max_cpu = force
        extradiv = 1
    else:
        max_cpu = None
        # without psutil, hyperthreaded CPU may overestimate physical count by factor of 2 (or more)
        if psutil is not None:
            max_cpu = psutil.cpu_count(logical=False)
            extradiv = 1
            if max_cpu is None:
                max_cpu = psutil.cpu_count()
                extradiv = 2
        if max_cpu is None:
            max_cpu = os.cpu_count()
            extradiv = 2

    return max_cpu // extradiv


def get_mpi_count(path: Pathlike, force: int = None) -> int:

    max_cpu = get_cpu_count(force)

    size = get_simsize(path)

    mpi_count = 1
    if size[2] == 1:  # 2D sim
        for i in range(max_cpu, 2, -1):
            mpi_count = max(math.gcd(size[1] // 2, i), mpi_count)
            if i < mpi_count:
                break
    else:  # 3D sim
        for i in range(max_cpu, 2, -1):
            mpi_count = max(math.gcd(size[2] // 2, i), mpi_count)
            if i < mpi_count:
                break

    return max(mpi_count, 1)


def get_simsize(path: Pathlike) -> T.Tuple[int, ...]:

    path = Path(path).expanduser()
    if path.is_dir():
        for suffix in (".h5", ".nc", ".dat"):
            fn = path / ("simsize" + suffix)
            if fn.is_file():
                break
    else:
        fn = path
    if not fn.is_file():
        raise FileNotFoundError(path)

    if fn.suffix == '.h5':
        return get_simsize_h5(fn)
    elif fn.suffix == '.nc':
        raise ValueError('TODO: implement NetCDF4')
    else:
        return get_simsize_raw(fn)
