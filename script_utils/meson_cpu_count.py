#!/usr/bin/env python3
import os
import argparse
from pathlib import Path
import struct
import math

try:
    import psutil
except ImportError:
    psutil = None
    # pip install psutil will improve CPU utilization.


def get_mpi_count(fn: Path, force: int = None) -> int:
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

    size = get_simsize(fn)
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

    return max(mpi_count, 2)


def get_simsize(fn: Path) -> tuple:
    fn = Path(fn).resolve().expanduser()
    if fn.stat().st_size != 12:
        raise ValueError(f"{fn} is not expected 12 bytes long")
    return struct.unpack("III", fn.open("rb").read(12))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fn", help="simsize.dat to read")
    p.add_argument("-f", "--force", help="force CPU count", type=int)
    P = p.parse_args()

    mpi_count = get_mpi_count(P.fn, P.force)

    # need end='' or you'll have to .strip() in Meson
    print(mpi_count, end="")
