from pathlib import Path
import math

from .utils import get_cpu_count
from .base import get_simsize
from .config import read_config


def get_mpi_count(path: Path, max_cpu: int) -> int:
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

    if not max_cpu:
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
