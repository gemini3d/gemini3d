from pathlib import Path
import typing as T
import numpy as np
from datetime import datetime

from .raw import get_simsize as get_simsize_raw

try:
    from .hdf import (
        get_simsize as get_simsize_h5,
        write_grid as write_grid_h5,
        write_state as write_state_h5,
        write_Efield as write_Efield_h5,
        write_precip as write_precip_h5,
    )
except ModuleNotFoundError:
    get_simsize_h5 = write_grid_h5 = write_state_h5 = None
try:
    from .nc4 import (
        get_simsize as get_simsize_nc,
        write_grid as write_grid_nc,
        write_state as write_state_nc,
        write_Efield as write_Efield_nc,
        write_precip as write_precip_nc,
    )
except ModuleNotFoundError:
    get_simsize_nc = write_grid_nc = write_state_nc = None

Pathlike = T.Union[str, Path]


def get_simsize(path: Pathlike) -> T.Tuple[int, ...]:

    path = Path(path).expanduser().resolve()
    if path.is_dir():
        for suffix in (".h5", ".nc", ".dat"):
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
        if get_simsize_h5 is None:
            raise ModuleNotFoundError("pip install h5py")
        return get_simsize_h5(fn)
    elif fn.suffix == ".nc":
        if get_simsize_nc is None:
            raise ModuleNotFoundError("pip install netcdf4")
        return get_simsize_nc(fn)
    else:
        return get_simsize_raw(fn)


def write_grid(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    """ writes grid to disk

    Parameters
    ----------

    p: dict
        simulation parameters
    xg: dict
        grid values
    """

    if p["format"] in ("hdf5", "h5"):
        if write_grid_h5 is None:
            raise ImportError("pip install h5py")
        write_grid_h5(p, xg)
    elif p["format"] in ("netcdf", "nc"):
        if write_grid_nc is None:
            raise ImportError("pip install netcdf4")
        write_grid_nc(p, xg)
    else:
        raise ValueError(f'unknown file format {p["format"]}')


def write_Efield(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    """ writes grid to disk

    Parameters
    ----------

    p: dict
        simulation parameters
    xg: dict
        grid values
    """

    if p["format"] in ("hdf5", "h5"):
        if write_Efield_h5 is None:
            raise ImportError("pip install h5py")
        write_Efield_h5(p, xg)
    elif p["format"] in ("netcdf", "nc"):
        if write_Efield_nc is None:
            raise ImportError("pip install netcdf4")
        write_Efield_nc(p, xg)
    else:
        raise ValueError(f'unknown file format {p["format"]}')


def write_precip(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    """ writes grid to disk

    Parameters
    ----------

    p: dict
        simulation parameters
    xg: dict
        grid values
    """

    if p["format"] in ("hdf5", "h5"):
        if write_precip_h5 is None:
            raise ImportError("pip install h5py")
        write_precip_h5(xg)
    elif p["format"] in ("netcdf", "nc"):
        if write_precip_nc is None:
            raise ImportError("pip install netcdf4")
        write_precip_nc(xg)
    else:
        raise ValueError(f'unknown file format {p["format"]}')


def write_state(time: datetime, ns: np.ndarray, vs: np.ndarray, Ts: np.ndarray, out_dir: Path, file_format: str):
    """
     WRITE STATE VARIABLE DATA.
    NOTE THAT WE don't write ANY OF THE ELECTRODYNAMIC
    VARIABLES SINCE THEY ARE NOT NEEDED TO START THINGS
    UP IN THE FORTRAN CODE.

    INPUT ARRAYS SHOULD BE TRIMMED TO THE CORRECT SIZE
    I.E. THEY SHOULD NOT INCLUDE GHOST CELLS
    """

    if file_format in ("hdf5", "h5"):
        if write_grid_h5 is None:
            raise ImportError("pip install h5py")
        write_state_h5(time, ns, vs, Ts, out_dir)
    elif file_format in ("netcdf", "nc"):
        if write_grid_nc is None:
            raise ImportError("pip install netcdf4")
        write_state_nc(time, ns, vs, Ts, out_dir)
    else:
        raise ValueError(f"unknown file format {file_format}")
