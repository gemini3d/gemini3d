#!/usr/bin/env python3
"""
python test_diffusion.py ../../../build/src/numerical/diffusion/test_diffusion1d_dirichlet.h5 -p

python test_diffusion.py ../../../build/src/numerical/diffusion/test_diffusion1d_neumann.h5 -p
"""
from pathlib import Path
import numpy as np
import h5py
import argparse
import sys


def read_diffusion1D(fn: Path, doplot: bool = False):
    """
    arrays are Ntime x Nx1

    Parameters
    ----------
    fn : Path
        HDF5 file to load
    doplot : bool, optional
        make plots
    """

    fn = Path(fn).expanduser().resolve()
    if not fn.is_file():
        print(fn, "not found", file=sys.stderr)
        raise SystemExit(77)

    with h5py.File(fn, "r") as f:
        lx1 = f["/lx1"][()]
        x1 = f["/x1"][:]
        t = f["/t"][:]
        TsEuler = f["/Euler/Ts"][:]
        TsBDF2 = f["/BDF2/Ts"][:]
        TsTrue = f["/True/Ts"][:]
        dirichTop = f["/flagdirichBottom"][()] == 1
        dirichBottom = f["/flagdirichTop"][()] == 1

    assert lx1 == TsEuler.shape[1] == TsBDF2.shape[1] == TsTrue.shape[1]

    if dirichBottom and dirichTop:
        assert np.isclose(
            TsEuler[12, 12], 0.770938954253086
        ), f"1-D Euler diffusion:dirichlet accuracy: expected 0.770939, got {TsEuler[12, 12]:.6f}"
        assert np.isclose(TsBDF2[12, 12], 0.763236513549944), "1-D BDF2 diffusion:dirichlet accuracy"
        assert np.isclose(TsTrue[12, 12], 0.763014494788105), "1-D true diffusion:dirichlet accuracy"
        print("OK: diffusion1d:dirichlet BCS")
    elif not dirichBottom and not dirichTop:
        # Intel 19.1, GCC 7.5:  got 0.966898
        assert np.isclose(
            TsEuler[12, 12], 0.966897705905538
        ), f"1-D Euler diffusion:neumann accuracy: expected 0.966897, got {TsEuler[12, 12]:.6f}"
        assert np.isclose(
            TsBDF2[12, 12], 0.9545709288834086
        ), f"1-D BDF2 diffusion:neumann accuracy: expected 0.954571, got {TsBDF2[12, 12]:.6f}"

    if not doplot:
        return

    fg = figure(figsize=(12, 5))
    ax = fg.subplots(1, 3, sharey=True, sharex=True)

    h = ax[0].pcolormesh(x1[2:-2], t, TsEuler, cmap="bwr")
    fg.colorbar(h, ax=ax[0])
    ax[0].set_xlabel("time [sec]")
    ax[0].set_ylabel("displacement [m]")
    ax[0].set_title("1-D diffusion (Backward Euler)")

    h = ax[1].pcolormesh(x1[2:-2], t, TsBDF2, cmap="bwr")
    fg.colorbar(h, ax=ax[1])
    ax[1].set_title("1-D diffusion (TRBDF2)")

    h = ax[2].pcolormesh(x1[2:-2], t, TsTrue, cmap="bwr")
    fg.colorbar(h, ax=ax[2])
    ax[2].set_title("1-D diffusion (analytic)")

    fg.suptitle(fn.name)


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
