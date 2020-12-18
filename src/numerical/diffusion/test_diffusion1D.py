#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import h5py
import argparse
import sys


def read_diffusion1D(fn: Path, doplot: bool = False):

    fn = Path(fn).expanduser()
    if not fn.is_file():
        print(fn, "not found", file=sys.stderr)
        raise SystemExit(77)

    with h5py.File(fn, "r") as f:
        lt = f["/lt"][()]
        lx1 = f["/lx1"][()]
        x1 = f["/x1"][:]
        t = np.empty(lt)
        TsEuler = np.empty((lx1, lt))
        TsBDF2 = np.empty((lx1, lt))
        Tstrue = np.empty((lx1, lt))
        for i in range(lt):
            t[i] = f[f"/t{i+1:04d}"][()]
            TsEuler[:, i] = f[f"/TsEuler{i+1:04d}"][:]
            TsBDF2[:, i] = f[f"/TsBDF2{i+1:04d}"][:]
            Tstrue[:, i] = f[f"/TsTrue{i+1:04d}"][:]

    assert np.isclose(
        TsEuler[12, 12], 0.770938954253086, rtol=0.01
    ), f"1-D Euler diffusion accuracy: truth: {TsEuler[12, 12]}"
    assert np.isclose(TsBDF2[12, 12], 0.763236513549944, rtol=0.01), "1-D BDF2 diffusion accuracy"
    assert np.isclose(Tstrue[12, 12], 0.763014494788105, rtol=1e-4), "1-D true diffusion accuracy"

    print("OK: 1d diffusion")

    if not doplot:
        return

    fg = figure()
    ax = fg.subplots(1, 3, sharey=True, sharex=True)

    h = ax[0].pcolormesh(t, x1[2:-2], TsEuler, cmap="bwr")
    fg.colorbar(h, ax=ax[0])
    ax[0].set_xlabel("time [sec]")
    ax[0].set_ylabel("displacement [m]")
    ax[0].set_title("1-D diffusion (Backward Euler)")

    h = ax[1].pcolormesh(t, x1[2:-2], TsBDF2, cmap="bwr")
    fg.colorbar(h, ax=ax[1])
    ax[1].set_title("1-D diffusion (TRBDF2)")

    h = ax[2].pcolormesh(t, x1[2:-2], Tstrue, cmap="bwr")
    fg.colorbar(h, ax=ax[2])
    ax[2].set_title("1-D diffusion (analytic)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("file")
    p.add_argument("-p", "--plot", help="make plots", action="store_true")
    P = p.parse_args()

    if P.plot:
        from matplotlib.pyplot import figure, show

    read_diffusion1D(P.file, P.plot)

    if P.plot:
        show()

    # with fn.open("r") as f:
    #     lt, lx1 = np.fromfile(f, int, 2, sep=" ")
    #     x1 = np.fromfile(f, float, lx1 + 4, sep=" ")
    #
    #     for i in range(lt):
    #         t[i] = np.fromfile(f, float, 1, sep=" ")
    #         TsEuler[:, i] = np.fromfile(f, float, lx1, sep=" ")
    #         TsBDF2[:, i] = np.fromfile(f, float, lx1, sep=" ")
    #         Tstrue[:, i] = np.fromfile(f, float, lx1, sep=" ")
