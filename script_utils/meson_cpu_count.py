#!/usr/bin/env python
import os
import argparse
from pathlib import Path
import numpy as np
import math
import typing


def cpu_count(min_cpu: int = 1) -> int:
    """
    This has the usual imperfection of reporting logical cores as well as physical
    cores. It may also report CPUs not actually available for your PID.

    We allow min_cpu to be larger than the actual CPU count for cases where a single-core CPU exists.
    For example, CI but we need to do tests that fail without a minimum image count.
    This scenario is inefficient but is usually just for CI.
    """
    if not min_cpu or min_cpu < 1:
        min_cpu = 1
    else:
        min_cpu = int(min_cpu)

    max_cpu = os.cpu_count()
    if not max_cpu:
        max_cpu = min_cpu
    else:
        max_cpu = max(max_cpu, min_cpu)

    return max_cpu


def get_simsize(fn: Path) -> typing.List[int]:
    fn = Path(fn).expanduser().resolve(strict=True)
    with fn.open("rb") as f:
        return np.fromfile(f, np.int32, 3).tolist()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fn", help="simsize.dat to read")
    p.add_argument(
        "min_cpu",
        help="minimum number of CPU cores to report. Some applications like MPI need at least some value to work.",
        type=int,
        nargs="?",
        default=1,
    )
    P = p.parse_args()

    max_cpu = cpu_count(P.min_cpu)

    size = get_simsize(P.fn)
    # print(size)
    if size[2] == 1:  # 2D sim
        mpi_count = math.gcd(size[1]//2, max_cpu)
    else:  # 3D sim
        mpi_count = math.gcd(size[2]//2, max_cpu)

    # need end='' or you'll have to .strip() in Meson
    print(mpi_count, end='')
