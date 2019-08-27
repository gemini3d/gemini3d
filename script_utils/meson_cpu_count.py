#!/usr/bin/env python3
import os
import argparse
from pathlib import Path
import numpy as np
import math
import typing

try:
    import psutil
except ImportError:
    psutil = None
    # pip install psutil will improve CPU utilization.


def get_simsize(fn: Path) -> typing.List[int]:
    fn = Path(fn).expanduser().resolve(strict=True)
    with fn.open("rb") as f:
        return np.fromfile(f, np.int32, 3).tolist()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fn", help="simsize.dat to read")
    p.add_argument("-f", "--force", help="force CPU count", type=int)
    P = p.parse_args()

    if P.force:
        max_cpu = P.force
        extradiv = 1
    else:
        # without psutil, hyperthreaded CPU may overestimate physical count by factor of 2 (or more)
        if psutil is not None:
            max_cpu = psutil.cpu_count(logical=False)
            extradiv = 1
        else:
            max_cpu = os.cpu_count()
            extradiv = 2

    size = get_simsize(P.fn)
    # need at least 2 images for MPI to function for Gemini
    mpi_count = 2
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

    mpi_count //= extradiv

    mpi_count = max(mpi_count, 2)

    # need end='' or you'll have to .strip() in Meson
    print(mpi_count, end="")
