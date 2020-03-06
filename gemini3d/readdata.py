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
from .base import get_simsize

try:
    from . import hdf
except ModuleNotFoundError:
    hdf = None

try:
    from . import nc4
except ModuleNotFoundError:
    nc4 = None

FILE_FORMATS = (".h5", ".nc", ".dat")


def readgrid(path: Path) -> typing.Dict[str, np.ndarray]:

    path = Path(path).expanduser().resolve()

    fn = get_grid_filename(path)

    if fn.suffix == ".dat":
        grid = raw.readgrid(fn)
    elif fn.suffix == ".h5":
        if hdf is None:
            raise ModuleNotFoundError("pip install h5py")
        grid = hdf.readgrid(fn)
    elif fn.suffix == ".nc":
        if nc4 is None:
            raise ModuleNotFoundError("pip install netcdf4")
        grid = nc4.readgrid(fn)
    else:
        raise ValueError(f"Unknown file type {fn}")

    return grid


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

    wavelength = [
        "3371",
        "4278",
        "5200",
        "5577",
        "6300",
        "7320",
        "10400",
        "3466",
        "7774",
        "8446",
        "3726",
        "LBH",
        "1356",
        "1493",
        "1304",
    ]

    fn = Path(fn).expanduser()
    fn_aurora = fn.parent / "aurmaps" / fn.name
    fn_Efield = fn.parent / "Efield_inputs" / fn.name

    input_dir = fn.parent / "inputs"
    P = read_config(input_dir)
    P["lxs"] = get_simsize(input_dir)

    if fn.suffix == ".dat":
        if P["flagoutput"] == 1:
            dat = raw.loadframe3d_curv(fn, P["lxs"])
        elif P["flagoutput"] == 2:
            dat = raw.loadframe3d_curvavg(fn, P["lxs"])
        else:
            raise ValueError("TODO: need to handle this case, file a bug report.")

        if fn_aurora.is_file():
            dat.update(raw.loadglow_aurmap(fn_aurora, P["lxs"], len(wavelength)))
            dat["wavelength"] = wavelength

        if fn_Efield.is_file():
            dat.update(read_Efield(fn_Efield))
    elif fn.suffix == ".h5":
        if hdf is None:
            raise ModuleNotFoundError("pip install h5py")

        if P["flagoutput"] == 1:
            dat = hdf.loadframe3d_curv(fn, P["lxs"])
        elif P["flagoutput"] == 2:
            dat = hdf.loadframe3d_curvavg(fn, P["lxs"])
        else:
            raise ValueError("TODO: need to handle this case, file a bug report.")

        if fn_aurora.is_file():
            dat.update(hdf.loadglow_aurmap(fn_aurora))
            dat["wavelength"] = wavelength
    elif fn.suffix == ".nc":
        raise NotImplementedError("TODO: NetCDF4")
    else:
        raise ValueError(f"Unknown file type {fn}")

    return dat


def read_Efield(fn: Path) -> typing.Dict[str, typing.Any]:
    """ load Efield data "Efield_inputs"

    Parameters
    ----------
    fn: pathlib.Path
        filename for this timestep

    Returns
    -------
    dat: dict of np.ndarray
        electric field
    """

    fn = Path(fn).expanduser().resolve(strict=True)

    if fn.suffix == ".dat":
        E = raw.load_Efield(fn)
    elif fn.suffix == ".h5":
        if hdf is None:
            raise ModuleNotFoundError("pip install h5py")
        E = hdf.load_Efield(fn)
    elif fn.suffix == ".nc":
        raise NotImplementedError("TODO: NetCDF4")
    else:
        raise ValueError(f"Unknown file type {fn}")

    return E


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

    return readdata(get_frame_filename(simdir, time))


def get_frame_filename(simdir: Path, time: datetime) -> Path:
    """
    the frame filenames can have different file formats
    """

    simdir = Path(simdir).expanduser().resolve(True)

    stem = time.strftime("%Y%m%d") + f"_{time.hour*3600 + time.minute*60 + time.second:05d}." + f"{time.microsecond:06d}"[:5]

    for ext in FILE_FORMATS:
        for tick in ("0", "1"):
            fn = simdir / (stem + tick + ext)
            if fn.is_file():
                return fn

    raise FileNotFoundError(f"could not find data file in {simdir} at {time}")


def get_grid_filename(path: Path) -> Path:
    """ given a path or filename, return the full path to simgrid file """

    path = Path(path).expanduser().resolve()

    if path.is_dir():
        for p in (path, path / "inputs"):
            for suff in FILE_FORMATS:
                file = p / ("simgrid" + suff)
                if file.is_file():
                    return file
    elif path.is_file():
        name = path.name
        path = path.parent
        for p in (path, path / "inputs"):
            file = p / name
            if file.is_file():
                return file

    raise FileNotFoundError(f"could not find grid file in {path}")
