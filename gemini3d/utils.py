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

from .config import read_config
from .base import get_simsize

git = shutil.which("git")

Pathlike = T.Union[str, Path]

__all__ = ["gitrev", "get_cpu_count", "get_mpi_count"]


def gitrev() -> str:
    if not git:
        return ""

    return subprocess.check_output([git, "rev-parse", "--short", "HEAD"], universal_newlines=True).strip()


def get_cpu_count() -> int:
    """ get a physical CPU count

    Returns
    -------
    count: int
        detect number of physical CPU
    """

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


def get_mpi_count(path: Path) -> int:
    """ get appropriate MPI image count for problem shape

    Parameters
    ----------
    path: pathlib.Path
        simsize file

    Returns
    -------
    count: int
        detect number of physical CPU
    """
    path = Path(path).expanduser()

    max_cpu = get_cpu_count()

    # %% config.nml file or directory or simsize.h5?
    if path.is_dir():
        size = get_simsize(path)
    elif path.is_file():
        if path.suffix in (".h5", ".nc", ".dat"):
            size = get_simsize(path)
        elif path.suffix in (".ini", ".nml"):
            params = read_config(path)
            # OK to use indat_size because we're going to run a sim on this machine
            size = get_simsize(params["indat_size"])
    else:
        raise FileNotFoundError(f"{path} is not a file or directory")

    mpi_count = 1
    if size[2] == 1:
        # 2D sim
        for i in range(max_cpu, 2, -1):
            mpi_count = max(math.gcd(size[1], i), mpi_count)
            if i < mpi_count:
                break
    else:
        # 3D sim
        for i in range(max_cpu, 2, -1):
            mpi_count = max(math.gcd(size[2], i), mpi_count)
            if i < mpi_count:
                break

    return mpi_count
