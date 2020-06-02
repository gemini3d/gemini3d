"""
struct manpage:
https://docs.python.org/3/library/struct.html#struct-format-strings
"""
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import typing as T

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

FILE_FORMATS = [".h5", ".nc", ".dat"]


def readgrid(path: Path, file_format: str = None) -> T.Dict[str, np.ndarray]:

    fn = get_grid_filename(path)

    if not file_format:
        file_format = fn.suffix[1:]

    if file_format == "dat":
        grid = raw.readgrid(fn.with_suffix(".dat"))
    elif file_format == "h5":
        if hdf is None:
            raise ModuleNotFoundError("pip install h5py")
        grid = hdf.readgrid(fn.with_suffix(".h5"))
    elif file_format == "nc":
        if nc4 is None:
            raise ModuleNotFoundError("pip install netcdf4")
        grid = nc4.readgrid(fn.with_suffix(".nc"))
    else:
        raise ValueError(f"Unknown file type {fn}")

    return grid


def readdata(fn: Path, file_format: str = None, *, E0dir: Path = None) -> T.Dict[str, T.Any]:
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
    if E0dir:
        fn_Efield = E0dir / fn.name

    input_dir = fn.parent / "inputs"
    P = read_config(input_dir)
    P["lxs"] = get_simsize(input_dir)

    if not file_format:
        file_format = fn.suffix[1:]

    if file_format == "dat":
        if P["flagoutput"] == 1:
            dat = raw.loadframe3d_curv(fn, P["lxs"])
        elif P["flagoutput"] == 2:
            dat = raw.loadframe3d_curvavg(fn, P["lxs"])
        else:
            raise ValueError("TODO: need to handle this case, file a bug report.")

        if fn_aurora.is_file():
            dat.update(raw.loadglow_aurmap(fn_aurora, P["lxs"], len(wavelength)))
            dat["wavelength"] = wavelength

        if E0dir and fn_Efield.is_file():
            dat.update(read_Efield(fn_Efield))
    elif file_format == "h5":
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
    elif file_format == "nc":
        if nc4 is None:
            raise ModuleNotFoundError("pip install netcdf4")

        if P["flagoutput"] == 1:
            dat = nc4.loadframe3d_curv(fn, P["lxs"])
        elif P["flagoutput"] == 2:
            dat = nc4.loadframe3d_curvavg(fn, P["lxs"])
        else:
            raise ValueError("TODO: need to handle this case, file a bug report.")

        if fn_aurora.is_file():
            dat.update(nc4.loadglow_aurmap(fn_aurora))
            dat["wavelength"] = wavelength
    else:
        raise ValueError(f"Unknown file type {fn}")

    return dat


def read_Efield(fn: Path, file_format: str = None) -> T.Dict[str, T.Any]:
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

    if not file_format:
        file_format = fn.suffix[1:]

    if file_format == "h5":
        if hdf is None:
            raise ModuleNotFoundError("pip install h5py")
        E = hdf.read_Efield(fn)
    elif file_format == "nc":
        if nc4 is None:
            raise ModuleNotFoundError("pip install netcdf4")
        E = nc4.read_Efield(fn)
    elif file_format == "dat":
        E = raw.read_Efield(fn)
    else:
        raise ValueError(f"Unknown file type {fn}")

    return E


def read_precip(fn: Path, file_format: str = None) -> T.Dict[str, T.Any]:
    """ load precipitation to disk

    Parameters
    ----------
    fn: pathlib.Path
        path to precipitation file
    file_format: str
        file format to read

    Returns
    -------
    dat: dict
        precipitation
    """

    fn = Path(fn).expanduser().resolve(strict=True)

    if not file_format:
        file_format = fn.suffix[1:]

    if file_format == "h5":
        if hdf is None:
            raise ImportError("pip install h5py")
        dat = hdf.read_precip(fn)
    elif file_format == "nc":
        if nc4 is None:
            raise ImportError("pip install netcdf4")
        dat = nc4.read_precip(fn)
    else:
        raise ValueError(f"unknown file format {file_format}")

    return dat


def read_state(file: Path,) -> T.Dict[str, T.Any]:
    """
    load inital condition data
    """

    if file.suffix == ".h5":
        if hdf is None:
            raise ImportError("pip install h5py")
        dat = hdf.read_state(file)
    elif file.suffix == ".nc":
        if nc4 is None:
            raise ImportError("pip install netcdf4")
        dat = nc4.read_state(file)
    else:
        raise ValueError(f"unknown file format {file.suffix}")

    return dat


def datetime_range(start: datetime, stop: datetime, step: timedelta) -> T.List[datetime]:

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


def loadframe(simdir: Path, time: datetime, file_format: str = None) -> T.Dict[str, T.Any]:
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
    file_format: str, optional
        "hdf5", "nc" for hdf5 or netcdf4 respectively

    Returns
    -------
    dat: dict
        simulation output for this time step
    """

    return readdata(get_frame_filename(simdir, time, file_format), file_format)


def get_frame_filename(simdir: Path, time: datetime, file_format: str = None) -> Path:
    """
    the frame filenames can have different file formats
    """

    simdir = Path(simdir).expanduser().resolve(True)

    stem = (
        time.strftime("%Y%m%d")
        + f"_{time.hour*3600 + time.minute*60 + time.second:05d}."
        + f"{time.microsecond:06d}"[:5]
    )

    suffixes = [f".{file_format}"] if file_format else FILE_FORMATS

    for ext in suffixes:
        for tick in ("0", "1"):
            fn = simdir / (stem + tick + ext)
            if fn.is_file():
                return fn

    raise FileNotFoundError(f"could not find data file in {simdir} at {time}")


def get_grid_filename(path: Path) -> Path:
    """ given a path or filename, return the full path to simgrid file
    we don't override FILE_FORMATS to allow outputs from a prior sim in a different
    file format to be used in this sim.
    """

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
