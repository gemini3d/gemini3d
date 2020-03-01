"""
NetCDF4 file IO
"""

from netCDF4 import Dataset
from pathlib import Path
import typing as T

LSP = 7


def get_simsize(path: Path) -> T.Tuple[int, ...]:
    """
    get simulation size

    TODO: needs to be validated
    """
    path = Path(path).expanduser().resolve()

    with Dataset(path, 'r') as f:
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
