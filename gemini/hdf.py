from pathlib import Path
import typing
import numpy as np
import logging
import h5py
from datetime import datetime, timedelta

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


def load_Efield(fn: Path) -> typing.Dict[str, typing.Any]:
    """
    load Efield_inputs files that contain input electric field in V/m
    """

    E: typing.Dict[str, np.ndarray] = {}

    sizefn = fn.parent / "simsize.h5"  # NOT the whole sim simsize.dat
    with h5py.File(sizefn, "r") as f:
        E["Nlon"] = f["Nlon"][()]
        E["Nlat"] = f["Nlat"][()]

    gridfn = fn.parent / "simgrid.h5"  # NOT the whole sim simgrid.dat
    with h5py.File(gridfn, "r") as f:
        E["mlon"] = f["mlon"][:]
        E["mlat"] = f["mlat"][:]

    with h5py.File(fn, "r") as f:
        E["flagdirich"] = f["flagdirich"]
        for p in ("Exit", "Eyit", "Vminx1it", "Vmaxx1it"):
            E[p] = [("x2", "x3"), f[p][:]]
        for p in ("Vminx2ist", "Vmaxx2ist"):
            E[p] = [("x2",), f[p][:]]
        for p in ("Vminx3ist", "Vmaxx3ist"):
            E[p] = [("x3",), f[p][:]]

    return E


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
        dat["time"] = ymdhourdec2datetime(f["time/ymd"][0], f["time/ymd"][1], f["time/ymd"][2], f["time/UThour"][()])

        dat["ne"] = (("x1", "x2", "x3"), f["/nsall"][LSP - 1, :, :, :].transpose())

        dat["v1"] = (
            ("x1", "x2", "x3"),
            (f["/nsall"][:6, :, :, :].transpose() * f["/vs1all"][:6, :, :, :].transpose()).sum(axis=3) / dat["ne"][1],
        )

        dat["Ti"] = (
            ("x1", "x2", "x3"),
            (f["/nsall"][:6, :, :, :].transpose() * f["/Tsall"][:6, :, :, :].transpose()).sum(axis=3) / dat["ne"][1],
        )
        dat["Te"] = (("x1", "x2", "x3"), f["/Tsall"][LSP - 1, :, :, :].transpose())

        dat["J1"] = (("x1", "x2", "x3"), f["/J1all"][:].transpose())
        dat["J2"] = (("x1", "x2", "x3"), f["/J2all"][:].transpose())
        dat["J3"] = (("x1", "x2", "x3"), f["/J3all"][:].transpose())

        dat["v2"] = (("x1", "x2", "x3"), f["/v2avgall"][:].transpose())
        dat["v3"] = (("x1", "x2", "x3"), f["/v3avgall"][:].transpose())

        dat["Phitop"] = (("x2", "x3"), f["/Phiall"][:])

    return dat


def loadframe3d_curvavg(fn: Path) -> typing.Dict[str, typing.Any]:
    """
    end users should normally use loadframe() instead

    Parameters
    ----------
    fn: pathlib.Path
        filename of this timestep of simulation output
    """
    #    grid = readgrid(fn.parent / "inputs/simgrid.h5")
    #    dat = xarray.Dataset(
    #        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    #    )
    dat: typing.Dict[str, typing.Any] = {}

    with h5py.File(fn, "r") as f:
        dat["time"] = ymdhourdec2datetime(f["time/ymd"][0], f["time/ymd"][1], f["time/ymd"][2], f["/time/UThour"][()])

        dat["ne"] = [("x1", "x2", "x3"), f["/neall"][:].transpose(2, 0, 1)]
        dat["v1"] = [("x1", "x2", "x3"), f["/v1avgall"][:].transpose(2, 0, 1)]
        dat["Ti"] = [("x1", "x2", "x3"), f["/Tavgall"][:].transpose(2, 0, 1)]
        dat["Te"] = [("x1", "x2", "x3"), f["/TEall"][:].transpose(2, 0, 1)]
        dat["J1"] = [("x1", "x2", "x3"), f["/J1all"][:].transpose(2, 0, 1)]
        dat["J2"] = [("x1", "x2", "x3"), f["/J2all"][:].transpose(2, 0, 1)]
        dat["J3"] = [("x1", "x2", "x3"), f["/J3all"][:].transpose(2, 0, 1)]
        dat["v2"] = [("x1", "x2", "x3"), f["/v2avgall"][:].transpose(2, 0, 1)]
        dat["v3"] = [("x1", "x2", "x3"), f["/v3avgall"][:].transpose(2, 0, 1)]
        dat["Phitop"] = [("x2", "x3"), f["/Phiall"][:]]

    return dat


def loadglow_aurmap(fn: Path) -> typing.Dict[str, typing.Any]:
    """
    read the auroral output from GLOW

    Parameters
    ----------
    fn: pathlib.Path
        filename of this timestep of simulation output
    """

    with h5py.File(fn, "r") as h:
        dat = {"rayleighs": [("wavelength", "x2", "x3"), h["/aurora/iverout"][:]]}

    return dat


def ymdhourdec2datetime(year: int, month: int, day: int, hourdec: float) -> datetime:
    """
    convert year,month,day + decimal hour HH.hhh to time
    """
    return datetime(year, month, day, int(hourdec), int((hourdec * 60) % 60)) + timedelta(seconds=(hourdec * 3600) % 60)
