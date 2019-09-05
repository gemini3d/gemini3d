from pathlib import Path
import typing
import numpy as np
import logging
import h5py
from dateutil.parser import parse

LSP = 7


def readgrid(fn: Path) -> typing.Dict[str, np.ndarray]:
    """
    get simulation dimensions

    Parameters
    ----------
    fn: pathlib.Path
        filepath to simgrid.h5

    Returns
    -------
    grid: dict
        grid parameters
    """

    grid: typing.Dict[str, typing.Any] = {}

    if not fn.is_file():
        logging.error(f"{fn} grid file is not present. Will try to load rest of data.")
        return grid

    with h5py.File(fn, "r") as f:
        grid["lx"] = f["lxs"][:]

        for key in f.keys():
            grid[key] = f[key][:]

    return grid


def loadframe3d_curv(fn: Path) -> typing.Dict[str, typing.Any]:
    """
    end users should normally use loadframe() instead
    """

    #    grid = readgrid(fn.parent / "inputs/simgrid.h5")
    #    dat = xarray.Dataset(
    #        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    #    )

    dat: typing.Dict[str, typing.Any] = {}

    with h5py.File(fn, "r") as f:
        dat["time"] = parse(f["time"][()])

        dat["ne"] = [("x1", "x2", "x3"), f["ns"][:, :, :, LSP - 1]]

        dat["v1"] = [("x1", "x2", "x3"), (f["ns"][:, :, :, :6] * f["vs1"][:, :, :, :6]).sum(axis=3) / f["ns"][:, :, :, LSP - 1]]

        dat["Ti"] = [("x1", "x2", "x3"), (f["ns"][:, :, :, :6] * f["Ts"][:, :, :, :6]).sum(axis=3) / f["ns"][:, :, :, LSP - 1]]
        dat["Te"] = [("x1", "x2", "x3"), f["Ts"][:, :, :, LSP - 1]]

        for p in ("J1", "J2", "J3", "v2", "v3"):
            dat[p] = [("x1", "x2", "x3"), f[p][:]]

        dat["Phitop"] = [("x2", "x3"), f["Phitop"][:]]

    return dat


def loadframe3d_curvavg(fn: Path) -> typing.Dict[str, typing.Any]:
    """
    end users should normally use loadframe() instead

    Parameters
    ----------
    path: pathlib.Path
        filename of this timestep of simulation output
    """
    #    grid = readgrid(fn.parent / "inputs/simgrid.h5")
    #    dat = xarray.Dataset(
    #        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    #    )
    dat: typing.Dict[str, typing.Any] = {}

    with h5py.File(fn, "r") as f:
        dat["time"] = parse(f["time"][()])

        for p in ("ne", "v1", "Ti", "Te", "J1", "J2", "J3", "v2", "v3"):
            dat[p] = [("x1", "x2", "x3"), f[p][:]]

        dat["Phitop"] = [("x2", "x3"), f["Phitop"][:]]

    return dat
