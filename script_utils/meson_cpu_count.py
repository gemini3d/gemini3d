#!/usr/bin/env python
import os
import argparse
from pathlib import Path
import numpy as np
import math
import typing


def get_simsize(fn: Path) -> typing.List[int]:
    fn = Path(fn).expanduser().resolve(strict=True)
    with fn.open("rb") as f:
        return np.fromfile(f, np.int32, 3).tolist()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fn", help="simsize.dat to read")
    P = p.parse_args()

    max_cpu = os.cpu_count()

    size = get_simsize(P.fn)
    # print(size)
    if size[2] == 1:  # 2D sim
        mpi_count = math.gcd(size[1] // 2, max_cpu)
    else:  # 3D sim
        mpi_count = math.gcd(size[2] // 2, max_cpu)

    # need at least 2 images for MPI to function for Gemini
    mpi_count = max(mpi_count, 2)

    # need end='' or you'll have to .strip() in Meson
    print(mpi_count, end="")
