#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import argparse
import sys


def test_diffusion1D(fn: Path, doplot: bool = False):

    fn = Path(fn).expanduser()
    if not fn.is_file():
        print(fn, "not found", file=sys.stderr)
        raise SystemExit(77)

    with fn.open("r") as f:
        lt, lx1 = np.fromfile(f, int, 2, sep=" ")

        x1 = np.fromfile(f, float, lx1 + 4, sep=" ")

        # %M=fscanf(fid,'%f',lx1*5)';
        # %M=reshape(M,[lx1 5])';
        # %b=fscanf(fid,'%f',lx1)';

        Ts = np.empty((lx1, lt))
        t = np.empty(lt)
        for i in range(lt):
            t[i] = np.fromfile(f, float, 1, sep=" ")
            Ts[:, i] = np.fromfile(f, float, lx1, sep=" ")

    assert np.isclose(Ts[12, -1], 0.2757552094055), "1-D diffusion accuracy"

    if doplot:
        fg = figure()
        ax = fg.gca()
        hi = ax.pcolormesh(t, x1[2:-2], Ts)
        fg.colorbar(hi, ax=ax)
        ax.set_xlabel("time [sec]")
        ax.set_ylabel("displacement [m]")
        ax.set_title("1-D diffusion (vs. time)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("file")
    p.add_argument("-p", "--plot", help="make plots", action="store_true")
    P = p.parse_args()

    if P.plot:
        from matplotlib.pyplot import figure, show

    test_diffusion1D(P.file, P.plot)

    if P.plot:
        show()
