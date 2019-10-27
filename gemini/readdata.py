"""
struct manpage:
https://docs.python.org/3/library/struct.html#struct-format-strings
"""
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import typing

from . import raw
from .config import read_config

try:
    from . import hdf
except ModuleNotFoundError:
    hdf = None


def readgrid(fn: Path) -> typing.Dict[str, np.ndarray]:

    if fn.is_dir():
        fn = fn / "inputs/simgrid.h5"
        if not fn.is_file():
            fn = fn.with_suffix(".dat")

    if fn.suffix == ".dat":
        return raw.readgrid(fn)
    else:
        if hdf is None:
            raise ModuleNotFoundError("pip install h5py")
        return hdf.readgrid(fn)


def readdata(fn: Path) -> typing.Dict[str, typing.Any]:
    """
    knowing the filename for a simulation time step, read the data for that time step

    Parameters
    ----------
    fn: pathlib.Path
        filename for this timestep

    Returns
    -------
    dat: dict
        simulation outputs as numpy.ndarray
    """
    fn = Path(fn).expanduser()
    P = read_config(fn.parent / "inputs")

    if fn.suffix == ".dat":
        if P["flagoutput"] == 1:
            dat = raw.loadframe3d_curv(fn, P["lxs"])
        elif P["flagoutput"] == 2:
            dat = raw.loadframe3d_curvavg(fn, P["lxs"])
        else:
            raise ValueError("TODO: need to handle this case, file a bug report.")
    else:
        if hdf is None:
            raise ModuleNotFoundError("pip install h5py")
        if P["flagoutput"] == 1:
            dat = hdf.loadframe3d_curv(fn)
        elif P["flagoutput"] == 2:
            dat = hdf.loadframe3d_curvavg(fn)
        else:
            raise ValueError("TODO: need to handle this case, file a bug report.")
    return dat


def datetime_range(start: datetime, stop: datetime, step: timedelta) -> typing.List[datetime]:

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


def loadframe(simdir: Path, time: datetime) -> typing.Dict[str, typing.Any]:
    """
    This is what users should normally use.
    load a frame of simulation data, automatically selecting the correct
    functions based on simulation parameters

    Parameters
    ----------
    simdir: pathlib.Path
        top-level directory of simulation output
    time: datetime.datetime
        time to load from simulation output

    Returns
    -------
    dat: dict
        simulation output for this time step
    """
    simdir = Path(simdir).expanduser().resolve(True)
    # %% datfn

    t = time
    stem = f"{t.year}{t.month:02d}{t.day:02d}_{t.hour*3600 + t.minute*60 + t.second:05d}.00000"
    for ext in ("0.h5", "1.h5", "0.dat", "1.dat"):
        datfn = simdir / (stem + ext)
        if datfn.is_file():
            break

    dat = readdata(datfn)

    return dat
