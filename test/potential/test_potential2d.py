#!/usr/bin/env python3
from pathlib import Path
import argparse

import numpy as np
import h5py


def read_potential2D(fn: str | Path) -> dict:

    fn = Path(fn).expanduser()

    with h5py.File(fn, "r") as f:
        v = {
            "lx1": f["/lx1"][()],
            "lx2": f["/lx2"][()],
            "lx3": f["/lx3"][()],
            "x1": f["/x1"][:],
            "x2": f["/x2"][:],
            "x3": f["/x3"][:],
            "Phi": f["/Phi"][:],
            "Phi2": f["/Phi2squeeze"][:],
            "Phitrue": f["/Phitrue"][:],
        }

    return v


def check_pot2d(v: dict) -> None:

    assert np.isclose(v["Phi2"][12, 12], 0.00032659, 1e-3), "Potential 2d accuracy"
    assert v["lx1"] == v["x1"].size
    assert v["lx2"] == v["x2"].size
    assert v["lx3"] == v["x3"].size


def plot_pot2d(v: dict) -> None:

    fg = figure(figsize=(15, 6))
    ax = fg.subplots(1, 3, sharey=True)
    h = ax[0].pcolormesh(v["x2"], v["x3"], v["Phi"])
    fg.colorbar(h, ax=ax[0])
    ax[0].set_ylabel("distance [m]")
    ax[0].set_xlabel("distance [m]")
    ax[0].set_title("2D potential (polarization)")

    h = ax[1].pcolormesh(v["x2"], v["x3"], v["Phi2"])
    fg.colorbar(h, ax=ax[1])
    ax[1].set_title("2D potential (static)")

    h = ax[2].pcolormesh(v["x2"], v["x3"], v["Phitrue"])
    fg.colorbar(h, ax=ax[2])
    ax[2].set_title("2D potential (analytical)")


def read_potential2D_old(fn: Path):
    """
    this isn't used anymore, but for reference shows hwo the raw binary files from Fortran used to be read.
    """
    with fn.open("r") as f:
        (lx2,) = np.fromfile(f, int, 1, sep=" ")
        x2 = np.fromfile(f, float, lx2, sep=" ")
        (lx3,) = np.fromfile(f, int, 1, sep=" ")
        x3 = np.fromfile(f, float, lx3, sep=" ")
        Phi = np.fromfile(f, float, lx2 * lx3, sep=" ").reshape((lx2, lx3))
        Phi2 = np.fromfile(f, float, lx2 * lx3, sep=" ").reshape((lx2, lx3))
        Phitrue = np.fromfile(f, float, lx2 * lx3, sep=" ").reshape((lx2, lx3))

    assert np.isclose(Phi2[12, 12], 0.000327, 1e-5), "Potential 2d accuracy"
    assert lx2 == x2.size
    assert lx3 == x3.size


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Test 2D potential solver output from Fortran MUMPS solver"
    )
    p.add_argument("file")
    p.add_argument("-p", "--plot", help="make plots", action="store_true")
    P = p.parse_args()

    v = read_potential2D(P.file)
    check_pot2d(v)

    if P.plot:
        from matplotlib.pyplot import figure, show

        plot_pot2d(v)
        show()
