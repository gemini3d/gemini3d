#!/usr/bin/env python3
import numpy as np
from pathlib import Path
import argparse
import sys


def read_potential2D(fn: Path, doplot: bool = False):

    fn = Path(fn).expanduser()
    if not fn.is_file():
        print(fn, "not found", file=sys.stderr)
        raise SystemExit(77)

    with fn.open("r") as f:
        (lx2,) = np.fromfile(f, int, 1, sep=" ")
        x2 = np.fromfile(f, float, lx2, sep=" ")
        (lx3,) = np.fromfile(f, int, 1, sep=" ")
        x3 = np.fromfile(f, float, lx3, sep=" ")
        Phi = np.fromfile(f, float, lx2 * lx3, sep=" ").reshape((lx2, lx3))
        Phi2 = np.fromfile(f, float, lx2 * lx3, sep=" ").reshape((lx2, lx3))
        Phitrue = np.fromfile(f, float, lx2 * lx3, sep=" ").reshape((lx2, lx3))
    assert np.isclose(Phi2[12, 12], 0.000327, 1e-5), "Potential 2d accuracy"

    if not doplot:
        return

    fg = figure(figsize=(15, 6))
    ax = fg.subplots(1, 3, sharey=True)
    h = ax[0].pcolormesh(x2, x3, Phi)
    fg.colorbar(h, ax=ax[0])
    ax[0].set_ylabel("distance [m]")
    ax[0].set_xlabel("distance [m]")
    ax[0].set_title("2D potential (polarization)")

    h = ax[1].pcolormesh(x2, x3, Phi2)
    fg.colorbar(h, ax=ax[1])
    ax[1].set_title("2D potential (static)")

    h = ax[2].pcolormesh(x2, x3, Phitrue)
    fg.colorbar(h, ax=ax[2])
    ax[2].set_title("2D potential (analytical)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("file")
    p.add_argument("-p", "--plot", help="make plots", action="store_true")
    P = p.parse_args()

    if P.plot:
        from matplotlib.pyplot import figure, show

    read_potential2D(P.file, P.plot)

    if P.plot:
        show()
