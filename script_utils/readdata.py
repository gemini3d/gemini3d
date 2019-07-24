import numpy as np
from pathlib import Path
from argparse import ArgumentParser
from datetime import datetime, timedelta
from typing import Tuple, Dict, Any

LSP = 7


def readsimsize(fn: Path):
    with fn.open("rb") as f:
        return np.fromfile(f, np.int32, 3)


def readdata(P: dict, fn: Path) -> Dict[str, np.ndarray]:

    dat: Dict[str, np.ndarray] = {}

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

        dat["Phitop"] = read2D(f, P["lxs"])

    if P["lxs"][1] == 1 or P["lxs"][2] == 1:
        dat["Ti"] = dat["Ti"].squeeze()
        dat["Te"] = dat["Te"].squeeze()

    return dat


def read4D(f, lsp: Tuple[int, int, int], lxs: int) -> np.ndarray:

    return np.fromfile(f, np.float64, np.prod(lxs) * lsp).reshape((lxs, lsp), order="F")


def read3D(f, lxs: Tuple[int, int, int]) -> np.ndarray:

    return np.fromfile(f, np.float64, np.prod(lxs)).reshape(*lxs, order="F")


def read2D(f, lxs: Tuple[int, int]) -> np.ndarray:

    return np.fromfile(f, np.float64, np.prod(lxs[1:])).reshape(*lxs[1:], order="F")


def readconfig(inifn: Path) -> Dict[str, Any]:
    inifn = Path(inifn).expanduser()

    P: Dict[str, Any] = {}

    with inifn.open("r") as f:
        date = list(map(int, f.readline().split()[0].split(",")))[::-1]
        sec = float(f.readline().split()[0])

        P["t0"] = datetime(*date) + timedelta(seconds=sec)  # type: ignore  # mypy bug

        P["tdur"] = float(f.readline().split()[0])

        P["dtout"] = float(f.readline().split()[0])

        P["f107a"], P["f107"], P["Ap"] = map(float, f.readline().split()[0].split(","))

        P["tcfl"] = float(f.readline().split()[0])
        P["Teinf"] = float(f.readline().split()[0])

        P["flagpot"] = int(f.readline().split()[0])
        P["flagperiodic"] = int(f.readline().split()[0])
        P["flagoutput"] = int(f.readline().split()[0])
        P["flagcap"] = int(f.readline().split()[0])

    return P


def loadframe(simdir: Path) -> Dict[str, np.ndarray]:
    simdir = Path(simdir).expanduser()

    P = readconfig(simdir / "inputs/config.ini")

    sizefn = simdir / "inputs/simsize.dat"
    P["lxs"] = readsimsize(sizefn)

    # %% datfn

    t = P["t0"].timetuple()
    datfn = simdir / f"{t[0]}{t[1]:02d}{t[2]:02d}_{t[3]*3600+t[4]*60+t[5]}.000001.dat"

    dat = readdata(P, datfn)

    return dat


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("simdir")
    p = p.parse_args()

    dat = loadframe(p.simdir)
