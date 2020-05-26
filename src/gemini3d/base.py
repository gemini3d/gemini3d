from pathlib import Path
import typing as T
import numpy as np
from datetime import datetime

from . import raw

try:
    from . import hdf
except ModuleNotFoundError:
    hdf = None

try:
    from . import nc4
except ModuleNotFoundError:
    nc4 = None

Pathlike = T.Union[str, Path]


def get_simsize(path: Pathlike) -> T.Tuple[int, ...]:
    """ get simulation dimenions """

    path = Path(path).expanduser().resolve()
    if path.is_dir():
        for suffix in [".h5", ".nc", ".dat"]:
            fn = path / ("simsize" + suffix)
            if fn.is_file():
                break
    else:
        fn = path
        if not fn.stem == "simsize":
            fn = path.parent / ("simsize" + path.suffix)
    if not fn.is_file():
        raise FileNotFoundError(path)

    if fn.suffix == ".h5":
        if hdf is None:
            raise ModuleNotFoundError("pip install h5py")
        return hdf.get_simsize(fn)
    elif fn.suffix == ".nc":
        if nc4 is None:
            raise ModuleNotFoundError("pip install netcdf4")
        return nc4.get_simsize(fn)
    else:
        return raw.get_simsize(fn)


def write_grid(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    """ writes grid to disk

    Parameters
    ----------

    p: dict
        simulation parameters
    xg: dict
        grid values

    NOTE: we use .with_suffix() in case file_format was overriden by user
    that allows writing NetCDF4 and HDF5 by scripts using same input files
    """

    p["indat_size"].parent.mkdir(parents=True, exist_ok=True)

    if p["format"] in ("hdf5", "h5"):
        if hdf is None:
            raise ImportError("pip install h5py")
        hdf.write_grid(p["indat_size"].with_suffix(".h5"), p["indat_grid"].with_suffix(".h5"), xg)
    elif p["format"] in ("netcdf", "nc"):
        if nc4 is None:
            raise ImportError("pip install netcdf4")
        nc4.write_grid(p["indat_size"].with_suffix(".nc"), p["indat_grid"].with_suffix(".nc"), xg)
    else:
        raise ValueError(f'unknown file format {p["format"]}')


def write_Efield(E: T.Dict[str, T.Any], outdir: Path, file_format: str):
    """ writes E-field to disk

    Parameters
    ----------

    E: dict
        E-field values
    outdir: pathlib.Path
        directory to write files into
    file_format: str
        requested file format to write
    """

    print("write E-field data to", outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if file_format in ("hdf5", "h5"):
        if hdf is None:
            raise ImportError("pip install h5py")
        hdf.write_Efield(outdir, E)
    elif file_format in ("netcdf", "nc"):
        if nc4 is None:
            raise ImportError("pip install netcdf4")
        nc4.write_Efield(outdir, E)
    else:
        raise ValueError(f"unknown file format {file_format}")


def write_precip(precip: T.Dict[str, T.Any], outdir: Path, file_format: str):
    """ writes precipitation to disk

    Parameters
    ----------
    precip: dict
        preicipitation values
    outdir: pathlib.Path
        directory to write files into
    file_format: str
        requested file format to write
    """

    print("write precipitation data to", outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if file_format in ("hdf5", "h5"):
        if hdf is None:
            raise ImportError("pip install h5py")
        hdf.write_precip(outdir, precip)
    elif file_format in ("netcdf", "nc"):
        if nc4 is None:
            raise ImportError("pip install netcdf4")
        nc4.write_precip(outdir, precip)
    else:
        raise ValueError(f"unknown file format {file_format}")


def write_state(
    time: datetime, ns: np.ndarray, vs: np.ndarray, Ts: np.ndarray, out_file: Path,
):
    """
     WRITE STATE VARIABLE DATA.
    NOTE THAT WE don't write ANY OF THE ELECTRODYNAMIC
    VARIABLES SINCE THEY ARE NOT NEEDED TO START THINGS
    UP IN THE FORTRAN CODE.

    INPUT ARRAYS SHOULD BE TRIMMED TO THE CORRECT SIZE
    I.E. THEY SHOULD NOT INCLUDE GHOST CELLS
    """

    if out_file.suffix == ".h5":
        if hdf is None:
            raise ImportError("pip install h5py")
        hdf.write_state(time, ns, vs, Ts, out_file.with_suffix(".h5"))
    elif out_file.suffix == ".nc":
        if nc4 is None:
            raise ImportError("pip install netcdf4")
        nc4.write_state(time, ns, vs, Ts, out_file.with_suffix(".nc"))
    else:
        raise ValueError(f"unknown file format {out_file.suffix}")
