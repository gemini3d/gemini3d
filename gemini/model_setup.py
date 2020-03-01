"""
setup a new simulation
"""
from pathlib import Path
import typing as T

from .config import read_nml
from .grid import makegrid_cart3d
from .plasma import equilibrium_state, equilibrium_resample, Efield_BCs2d, Efield_BCs3d, particles_BCs
from .hdf import write_state, write_grid

R = Path(__file__).parents[1]


def model_setup(path: Path, out_dir: Path):
    """
  top-level function to create a new simulation

  Parameters
  ----------

  path: pathlib.Path
      path (directory or full path) to config.nml
  out_dir: pathlib.Path
      directory to write simulation artifacts to
  """

    # %% read config.nml
    p = read_nml(path)
    p["out_dir"] = Path(out_dir).expanduser().resolve()

    # %% is this equilibrium or interpolated simulation
    if "eqdir" in p:
        model_setup_interp(p)
    else:
        model_setup_equilibrium(p)


def model_setup_equilibrium(p: T.Dict[str, T.Any]):
    # %% GRID GENERATION

    xg = makegrid_cart3d(p)

    write_grid(p, xg)

    # %% Equilibrium input generation
    [ns, Ts, vsx1] = equilibrium_state(p, xg)

    write_state(p["t0"], ns, vsx1, Ts, p["out_dir"], p["format"], p["realbits"])


def model_setup_interp(p: T.Dict[str, T.Any]):

    xg = makegrid_cart3d(p)

    equilibrium_resample(p, xg)

    # %% potential boundary conditions
    if "flagE0file" in p and p["flagE0file"]:
        if p["lxp"] == 1 or p["lyp"] == 1:
            Efield_BCs2d(p)
        else:  # 3D
            Efield_BCs3d(p)

    # %% aurora
    if "flagprecfile" in p and p["flagprecfile"]:
        particles_BCs(p)
