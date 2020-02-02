from pathlib import Path
import typing
import numpy as np
import logging
import struct
from datetime import datetime, timedelta
from functools import lru_cache

LSP = 7


@lru_cache()
def get_simsize(fn: Path) -> tuple:
    """
    get simulation dimensions from simsize.dat
    in the future, this would be in the .h5 HDF5 output.

    Parameters
    ----------
    fn: pathlib.Path
        filepath to simsize.dat

    Returns
    -------
    size: tuple of int, int, int
        3 integers telling simulation grid size
    """
    fn = Path(fn).expanduser()
    if fn.stat().st_size != 12:
        raise ValueError(f"{fn} is not expected 12 bytes long")
    return struct.unpack("III", fn.open("rb").read(12))


def readgrid(fn: Path) -> typing.Dict[str, np.ndarray]:
    """
    get simulation dimensions
    in the future, this would be in the .h5 HDF5 output.

    Parameters
    ----------
    fn: pathlib.Path
        filepath to simgrid.dat

    Returns
    -------
    grid: dict
        grid parameters
    """
    lxs = get_simsize(fn.parent / "simsize.dat")
    lgridghost = (lxs[0] + 4) * (lxs[1] + 4) * (lxs[2] + 4)
    gridsizeghost = [lxs[0] + 4, lxs[1] + 4, lxs[2] + 4]

    grid: typing.Dict[str, typing.Any] = {"lx": lxs}

    if not fn.is_file():
        logging.error(f"{fn} grid file is not present. Will try to load rest of data.")
        return grid

    with fn.open("r") as f:
        read = np.fromfile
        for i in (1, 2, 3):
            grid[f"x{i}"] = read(f, np.float64, lxs[i - 1] + 4)
            grid[f"x{i}i"] = read(f, np.float64, lxs[i - 1] + 1)
            grid[f"dx{i}b"] = read(f, np.float64, lxs[i - 1] + 3)
            grid[f"dx{i}h"] = read(f, np.float64, lxs[i - 1])
        for i in (1, 2, 3):
            grid[f"h{i}"] = read(f, np.float64, lgridghost).reshape(gridsizeghost)
        L = [lxs[0] + 1, lxs[1], lxs[2]]
        for i in (1, 2, 3):
            grid[f"h{i}x1i"] = read(f, np.float64, np.prod(L)).reshape(L)
        L = [lxs[0], lxs[1] + 1, lxs[2]]
        for i in (1, 2, 3):
            grid[f"h{i}x2i"] = read(f, np.float64, np.prod(L)).reshape(L)
        L = [lxs[0], lxs[1], lxs[2] + 1]
        for i in (1, 2, 3):
            grid[f"h{i}x3i"] = read(f, np.float64, np.prod(L)).reshape(L)
        for i in (1, 2, 3):
            grid[f"gx{i}"] = read(f, np.float64, np.prod(lxs)).reshape(lxs)
        for k in ("alt", "glat", "glon", "Bmag"):
            grid[k] = read(f, np.float64, np.prod(lxs)).reshape(lxs)
        grid["Bincl"] = read(f, np.float64, lxs[1] * lxs[2]).reshape(lxs[1:])
        grid["nullpts"] = read(f, np.float64, np.prod(lxs)).reshape(lxs)
        if f.tell() == fn.stat().st_size:  # not EOF
            return grid

        L = [lxs[0], lxs[1], lxs[2], 3]
        for i in (1, 2, 3):
            grid[f"e{i}"] = read(f, np.float64, np.prod(L)).reshape(L)
        for k in ("er", "etheta", "ephi"):
            grid[k] = read(f, np.float64, np.prod(L)).reshape(L)
        for k in ("r", "theta", "phi"):
            grid[k] = read(f, np.float64, np.prod(lxs)).reshape(lxs)
        if f.tell() == fn.stat().st_size:  # not EOF
            return grid

        for k in ("x", "y", "z"):
            grid[k] = read(f, np.float64, np.prod(lxs)).reshape(lxs)

    return grid


def load_Efield(fn: Path) -> typing.Dict[str, np.ndarray]:
    """
    load Efield_inputs files that contain input electric field in V/m
    """

    read = np.fromfile

    E: typing.Dict[str, np.ndarray] = {}

    sizefn = fn.parent / "simsize.dat"  # NOT the whole sim simsize.dat
    with sizefn.open("r") as f:
        E["Nlon"], E["Nlat"] = read(f, np.int32, 2)

    assert E['Nlon'] > 0, 'must have strictly positive number of longitude cells'
    assert E['Nlat'] > 0, 'must have strictly positive number of latitude cells'

    lxs = (0, E["Nlon"], E["Nlat"])

    gridfn = fn.parent / "simgrid.dat"  # NOT the whole sim simgrid.dat
    with gridfn.open("r") as f:
        E["mlon"] = read(f, np.float64, E["Nlon"])
        E["mlat"] = read(f, np.float64, E["Nlat"])

    assert ((E['mlat'] >= -90) & (E['mlat'] <= 90)).all(), f'impossible latitude, was file read correctly? {gridfn}'

    with fn.open("r") as f:
        """
        NOTE:
        this is mistakenly a float from Matlab
        to keep compatibility with old files, we left it as real64.
        New work should be using HDF5 instead of raw in any case.
        """
        E["flagdirich"] = read(f, np.float64, 1)
        for p in ("Exit", "Eyit", "Vminx1it", "Vmaxx1it"):
            E[p] = read2D(f, lxs)
        for p in ("Vminx2ist", "Vmaxx2ist"):
            E[p] = read(f, np.float64, E["Nlat"])
        for p in ("Vminx3ist", "Vmaxx3ist"):
            E[p] = read(f, np.float64, E["Nlon"])
        filesize = fn.stat().st_size
        if f.tell() != filesize:
            logging.error(f'{fn} size {filesize} != file read position {f.tell()}')

    return E


def loadframe3d_curv(fn: Path, lxs: typing.Sequence[int]) -> typing.Dict[str, typing.Any]:
    """
    end users should normally use loadframe() instead

    Parameters
    ----------
    fn: pathlib.Path
        filename of this timestep of simulation output
    lxs: list of int
        array dimension
    """

    #    grid = readgrid(fn.parent / "inputs/simgrid.dat")
    #    dat = xarray.Dataset(
    #        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    #    )

    dat: typing.Dict[str, typing.Any] = {}

    with fn.open("r") as f:
        dat["time"] = read_time(f)

        ns = read4D(f, LSP, lxs)
        dat["ne"] = (("x1", "x2", "x3"), ns[:, :, :, LSP - 1])

        vs1 = read4D(f, LSP, lxs)
        dat["v1"] = (("x1", "x2", "x3"), (ns[:, :, :, :6] * vs1[:, :, :, :6]).sum(axis=3) / dat["ne"][1])

        Ts = read4D(f, LSP, lxs)
        dat["Ti"] = (("x1", "x2", "x3"), (ns[:, :, :, :6] * Ts[:, :, :, :6]).sum(axis=3) / dat["ne"][1])
        dat["Te"] = (("x1", "x2", "x3"), Ts[:, :, :, LSP - 1].squeeze())

        for p in ("J1", "J2", "J3", "v2", "v3"):
            dat[p] = [("x1", "x2", "x3"), read3D(f, lxs)]

        dat["Phitop"] = [("x2", "x3"), read2D(f, lxs)]

    return dat


def loadframe3d_curvavg(fn: Path, lxs: typing.Sequence[int]) -> typing.Dict[str, typing.Any]:
    """
    end users should normally use loadframe() instead

    Parameters
    ----------
    fn: pathlib.Path
        filename of this timestep of simulation output
    lxs: list of int
        array dimension
    """
    #    grid = readgrid(fn.parent / "inputs/simgrid.dat")
    #    dat = xarray.Dataset(
    #        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    #    )
    dat: typing.Dict[str, typing.Any] = {}

    with fn.open("r") as f:
        dat["time"] = read_time(f)

        for p in ("ne", "v1", "Ti", "Te", "J1", "J2", "J3", "v2", "v3"):
            dat[p] = [("x1", "x2", "x3"), read3D(f, lxs)]

        dat["Phitop"] = [("x2", "x3"), read2D(f, lxs)]

    return dat


def read4D(f, lsp: int, lxs: typing.Sequence[int]) -> np.ndarray:
    """
    end users should normally use laodframe() instead
    """
    if not len(lxs) == 3:
        raise ValueError(f"lxs must have 3 elements, you have lxs={lxs}")

    return np.fromfile(f, np.float64, np.prod(lxs) * lsp).reshape((*lxs, lsp), order="F")


def read3D(f, lxs: typing.Sequence[int]) -> np.ndarray:
    """
    end users should normally use loadframe() instead
    """
    if not len(lxs) == 3:
        raise ValueError(f"lxs must have 3 elements, you have lxs={lxs}")

    return np.fromfile(f, np.float64, np.prod(lxs)).reshape(*lxs, order="F")


def read2D(f, lxs: typing.Sequence[int]) -> np.ndarray:
    """
    end users should normally use laodframe() instead
    """
    if not len(lxs) == 3:
        raise ValueError(f"lxs must have 3 elements, you have lxs={lxs}")

    return np.fromfile(f, np.float64, np.prod(lxs[1:])).reshape(*lxs[1:], order="F")


def loadglow_aurmap(f, lxs: typing.Sequence[int], lwave: int) -> typing.Dict[str, typing.Any]:
    """
    read the auroral output from GLOW
    """
    if not len(lxs) == 3:
        raise ValueError(f"lxs must have 3 elements, you have lxs={lxs}")
    raw = np.fromfile(f, np.float64, np.prod(lxs[1:]) * lwave).reshape(np.prod(lxs[1:]) * lwave, order="F")
    return {"rayleighs": [("wavelength", "x2", "x3"), raw]}


def read_time(f) -> datetime:
    t = np.fromfile(f, np.float64, 4)
    return datetime(int(t[0]), int(t[1]), int(t[2])) + timedelta(hours=t[3])
