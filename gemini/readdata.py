import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import typing

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


def readsimsize(fn: Path):
    with fn.open("rb") as f:
        return np.fromfile(f, np.int32, 3)


def readdata(P: typing.Dict[str, typing.Any], fn: Path) -> typing.Dict[str, np.ndarray]:

    dat: typing.Dict[str, np.ndarray] = {}

    with fn.open("rb") as f:
        t = np.fromfile(f, np.float64, 4)
        dat["t"] = datetime(int(t[0]), int(t[1]), int(t[2])) + timedelta(hours=t[3])

        dat["ne"] = read3D(f, P["lxs"])

        dat["v1"] = read3D(f, P["lxs"])

        dat["Ti"] = read3D(f, P["lxs"])
        dat["Te"] = read3D(f, P["lxs"])

        dat["J1"] = read3D(f, P["lxs"])
        dat["J2"] = read3D(f, P["lxs"])
        dat["J3"] = read3D(f, P["lxs"])

        dat["v2"] = read3D(f, P["lxs"])
        dat["v3"] = read3D(f, P["lxs"])

        dat["Phitop"] = read2D(f, P["lxs"])

    if P["lxs"][1] == 1 or P["lxs"][2] == 1:
        dat["Ti"] = dat["Ti"].squeeze()
        dat["Te"] = dat["Te"].squeeze()

    return dat


def read4D(f, lsp: typing.Tuple[int, int, int], lxs: int) -> np.ndarray:

    return np.fromfile(f, np.float64, np.prod(lxs) * lsp).reshape((lxs, lsp), order="F")


def read3D(f, lxs: typing.Tuple[int, int, int]) -> np.ndarray:

    return np.fromfile(f, np.float64, np.prod(lxs)).reshape(*lxs, order="F")


def read2D(f, lxs: typing.Tuple[int, int]) -> np.ndarray:

    return np.fromfile(f, np.float64, np.prod(lxs[1:])).reshape(*lxs[1:], order="F")


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


def loadframe(simdir: Path, time: datetime) -> typing.Dict[str, np.ndarray]:
    simdir = Path(simdir).expanduser()

    P = readconfig(simdir / "inputs/config.ini")

    sizefn = simdir / "inputs/simsize.dat"
    P["lxs"] = readsimsize(sizefn)

    # %% datfn

    t = time.timetuple()
    timename = f"{t[0]}{t[1]:02d}{t[2]:02d}_{t[3]*3600+t[4]*60+t[5]}.000000.dat"
    datfn = simdir / timename
    if not datfn.is_file():
        datfn = datfn.parent / (timename[:-10] + "000001.dat")

    dat = readdata(P, datfn)

    return dat
