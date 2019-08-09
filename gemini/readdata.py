"""
struct manpage:
https://docs.python.org/3/library/struct.html#struct-format-strings
"""
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import typing
from functools import lru_cache
import struct
import xarray

LSP = 7


def datetime_range(
    start: datetime, stop: datetime, step: timedelta
) -> typing.List[datetime]:

    """
    Generate range of datetime

    Parameters
    ----------
    start : datetime
        start time
    stop : datetime
        stop time
    step : timedelta
        time step

    Returns
    -------
    times : list of datetime
        times requested
    """
    return [start + i * step for i in range((stop - start) // step)]


@lru_cache()
def get_simsize(fn: Path) -> typing.Tuple[int, int, int]:
    return struct.unpack("III", Path(fn).expanduser().read_bytes())  # type: ignore


def readgrid(fn: Path) -> typing.Dict[str, np.ndarray]:
    lxs = get_simsize(fn.parent / "simsize.dat")
    lgridghost = (lxs[0] + 4) * (lxs[1] + 4) * (lxs[2] + 4)
    gridsizeghost = [lxs[0] + 4, lxs[1] + 4, lxs[2] + 4]

    grid = {}
    with fn.open("rb") as f:
        for i in (1, 2, 3):
            grid[f"x{i}"] = np.fromfile(f, np.float64, lxs[i - 1] + 4)
            grid[f"x{i}i"] = np.fromfile(f, np.float64, lxs[i - 1] + 1)
            grid[f"dx{i}b"] = np.fromfile(f, np.float64, lxs[i - 1] + 3)
            grid[f"dx{i}h"] = np.fromfile(f, np.float64, lxs[i - 1])
        for i in (1, 2, 3):
            grid[f"h{i}"] = np.fromfile(f, np.float64, lgridghost).reshape(
                gridsizeghost
            )

    return grid


def readdata(fn: Path) -> xarray.Dataset:

    P = readconfig(fn.parent / "inputs/config.ini")
    if P["flagoutput"] == 1:
        dat = loadframe3d_curv(fn)
    elif P["flagoutput"] == 2:
        dat = loadframe3d_curvavg(fn)
    else:
        raise NotImplementedError("TODO: need to handle this case, file a bug report.")
    return dat


def loadframe3d_curv(fn: Path) -> xarray.Dataset:
    P = readconfig(fn.parent / "inputs/config.ini")
    grid = readgrid(fn.parent / "inputs/simgrid.dat")
    dat = xarray.Dataset(
        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    )

    with fn.open("rb") as f:
        t = np.fromfile(f, np.float64, 4)
        dat.attrs["time"] = datetime(int(t[0]), int(t[1]), int(t[2])) + timedelta(
            hours=t[3]
        )

        ns = read4D(f, LSP, P["lxs"])
        dat["ne"] = (("x1", "x2", "x3"), ns[:, :, :, LSP - 1].squeeze())

        vs1 = read4D(f, LSP, P["lxs"])
        dat["v1"] = (
            ("x1", "x2", "x3"),
            (ns[:, :, :, :6] * vs1[:, :, :, :6]).sum(axis=3) / ns[:, :, :, LSP - 1],
        )

        Ts = read4D(f, LSP, P["lxs"])
        dat["Ti"] = (
            ("x1", "x2", "x3"),
            (ns[:, :, :, :6] * Ts[:, :, :, :6]).sum(axis=3) / ns[:, :, :, LSP - 1],
        )
        dat["Te"] = (("x1", "x2", "x3"), Ts[:, :, :, LSP - 1].squeeze())

        for p in ("J1", "J2", "J3", "v2", "v3"):
            dat[p] = (("x1", "x2", "x3"), read3D(f, P["lxs"]))

        dat["Phitop"] = (("x2", "x3"), read2D(f, P["lxs"]))

    return dat


def loadframe3d_curvavg(fn: Path) -> xarray.Dataset:
    P = readconfig(fn.parent / "inputs/config.ini")
    grid = readgrid(fn.parent / "inputs/simgrid.dat")
    dat = xarray.Dataset(
        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    )

    with fn.open("rb") as f:
        t = np.fromfile(f, np.float64, 4)
        dat.attrs["time"] = datetime(int(t[0]), int(t[1]), int(t[2])) + timedelta(
            hours=t[3]
        )

        for p in ("ne", "v1", "Ti", "Te", "J1", "J2", "J3", "v2", "v3"):
            dat[p] = (("x1", "x2", "x3"), read3D(f, P["lxs"]))

        dat["Phitop"] = (("x2", "x3"), read2D(f, P["lxs"]))

    if P["lxs"][1] == 1 or P["lxs"][2] == 1:
        dat["Ti"] = dat["Ti"].squeeze()
        dat["Te"] = dat["Te"].squeeze()

    return dat


def read4D(f, lsp: int, lxs: typing.Sequence[int]) -> np.ndarray:
    if not len(lxs) == 3:
        raise ValueError(f"lxs must have 3 elements, you have lxs={lxs}")

    return np.fromfile(f, np.float64, np.prod(lxs) * lsp).reshape(
        (*lxs, lsp), order="F"
    )


def read3D(f, lxs: typing.Sequence[int]) -> np.ndarray:
    if not len(lxs) == 3:
        raise ValueError(f"lxs must have 3 elements, you have lxs={lxs}")

    return np.fromfile(f, np.float64, np.prod(lxs)).reshape(*lxs, order="F")


def read2D(f, lxs: typing.Sequence[int]) -> np.ndarray:
    if not len(lxs) == 3:
        raise ValueError(f"lxs must have 3 elements, you have lxs={lxs}")

    return np.fromfile(f, np.float64, np.prod(lxs[1:])).reshape(*lxs[1:], order="F")


@lru_cache()
def readconfig(inifn: Path) -> typing.Dict[str, typing.Any]:
    inifn = Path(inifn).expanduser().resolve(strict=True)

    P: typing.Dict[str, typing.Any] = {}

    with inifn.open("r") as f:
        date = list(map(int, f.readline().split()[0].split(",")))[::-1]
        sec = float(f.readline().split()[0])

        P["t0"] = datetime(*date) + timedelta(seconds=sec)  # type: ignore  # mypy bug

        P["tdur"] = timedelta(seconds=float(f.readline().split()[0]))

        P["dtout"] = timedelta(seconds=float(f.readline().split()[0]))

        P["f107a"], P["f107"], P["Ap"] = map(float, f.readline().split()[0].split(","))

        P["tcfl"] = float(f.readline().split()[0])
        P["Teinf"] = float(f.readline().split()[0])

        P["flagpot"] = int(f.readline().split()[0])
        P["flagperiodic"] = int(f.readline().split()[0])
        P["flagoutput"] = int(f.readline().split()[0])
        P["flagcap"] = int(f.readline().split()[0])

    return P


def loadframe(simdir: Path, time: datetime) -> xarray.Dataset:
    simdir = Path(simdir).expanduser()

    P = readconfig(simdir / "inputs/config.ini")

    sizefn = simdir / "inputs/simsize.dat"
    P["lxs"] = get_simsize(sizefn)

    # %% datfn

    t = time.timetuple()
    timename = f"{t[0]}{t[1]:02d}{t[2]:02d}_{t[3]*3600+t[4]*60+t[5]}.000000.dat"
    datfn = simdir / timename
    if not datfn.is_file():
        datfn = datfn.parent / (timename[:-10] + "000001.dat")

    dat = readdata(datfn)

    return dat
