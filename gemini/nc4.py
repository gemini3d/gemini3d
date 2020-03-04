"""
NetCDF4 file IO
"""

from netCDF4 import Dataset
from pathlib import Path
import typing as T
import numpy as np
from datetime import datetime

LSP = 7


def get_simsize(path: Path) -> T.Tuple[int, ...]:
    """
    get simulation size

    TODO: needs to be validated
    """
    path = Path(path).expanduser().resolve()

    with Dataset(path, "r") as f:
        if "lxs" in f:
            lxs = f["lxs"][:]
        elif "lx1" in f:
            if f["lx1"].ndim > 0:
                lxs = (f["lx1"][:].squeeze()[()], f["lx2"][:].squeeze()[()], f["lx3"][:].squeeze()[()])
            else:
                lxs = (f["lx1"][()], f["lx2"][()], f["lx3"][()])
        else:
            raise KeyError(f"could not find '/lxs' or '/lx1' in {path.as_posix()}")

    return lxs


def readgrid(fn: Path) -> T.Dict[str, np.ndarray]:
    raise NotImplementedError("TODO: NetCDF4 simgrid.nc")


def write_grid(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    """ writes grid to disk

    Parameters
    ----------

    p: dict
        simulation parameters
    xg: dict
        grid values
    """
    raise NotImplementedError("TODO: NetCDF4")


def write_state(time: datetime, ns: np.ndarray, vs: np.ndarray, Ts: np.ndarray, out_dir: Path):
    """
     WRITE STATE VARIABLE DATA.
    NOTE THAT WE don't write ANY OF THE ELECTRODYNAMIC
    VARIABLES SINCE THEY ARE NOT NEEDED TO START THINGS
    UP IN THE FORTRAN CODE.

    INPUT ARRAYS SHOULD BE TRIMMED TO THE CORRECT SIZE
    I.E. THEY SHOULD NOT INCLUDE GHOST CELLS
    """

    fn = out_dir / "initial_conditions.h5"
    print("write", fn)

    with Dataset(fn, "w") as f:
        f["/ymd"] = [time.year, time.month, time.day]
        f["/UTsec"] = time.hour * 3600 + time.minute * 60 + time.second + time.microsecond / 1e6
        f["/ns"] = ns
        f["/vsx1"] = vs
        f["/Ts"] = Ts


def write_Efield(p: T.Dict[str, T.Any], E: T.Dict[str, np.ndarray]):
    """
    write Efield to disk
    """
    raise NotImplementedError("TODO: NetCDF")


def write_precip(E: T.Dict[str, np.ndarray]):
    """
    write precipitation to disk
    """
    raise NotImplementedError("TODO: NetCDF")
