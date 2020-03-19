"""
NetCDF4 file IO
"""

from netCDF4 import Dataset
from pathlib import Path
import typing as T
import numpy as np
from datetime import datetime

from .utils import datetime2ymd_hourdec

LSP = 7


def get_simsize(path: Path) -> T.Tuple[int, ...]:
    """
    get simulation size
    """

    path = Path(path).expanduser().resolve()

    with Dataset(path, "r") as f:
        if "lxs" in f.variables:
            lxs = f["lxs"][:]
        elif "lx" in f.variables:
            lxs = f["lx"][:]
        elif "lx1" in f.variables:
            if f["lx1"].ndim > 0:
                lxs = (f["lx1"][:].squeeze()[()], f["lx2"][:].squeeze()[()], f["lx3"][:].squeeze()[()])
            else:
                lxs = (f["lx1"][()], f["lx2"][()], f["lx3"][()])
        else:
            raise KeyError(f"could not find 'lxs', 'lx' or 'lx1' in {path.as_posix()}")

    return lxs.data


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
    (p["out_dir"] / "inputs").mkdir(parents=True, exist_ok=True)

    fn = p["out_dir"] / "inputs/simsize.nc"
    print("write_grid:", fn)
    with Dataset(fn, "w") as f:
        f.createDimension("length", len(xg["lx"]))
        g = f.createVariable("lx", np.int32, ("length",))
        g[:] = xg["lx"]

    fn = p["out_dir"] / "inputs/simgrid.nc"
    print("write_grid:", fn)
    Ng = 4  # number of ghost cells

    with Dataset(fn, "w") as f:
        f.createDimension("x1ghost", xg["lx"][0] + Ng)
        f.createDimension("x1d", xg["lx"][0] + Ng - 1)
        f.createDimension("x1i", xg["lx"][0] + 1)
        f.createDimension("x1", xg["lx"][0])

        f.createDimension("x2ghost", xg["lx"][1] + Ng)
        f.createDimension("x2d", xg["lx"][1] + Ng - 1)
        f.createDimension("x2i", xg["lx"][1] + 1)
        f.createDimension("x2", xg["lx"][1])

        f.createDimension("x3ghost", xg["lx"][2] + Ng)
        f.createDimension("x3d", xg["lx"][2] + Ng - 1)
        f.createDimension("x3i", xg["lx"][2] + 1)
        f.createDimension("x3", xg["lx"][2])

        f.createDimension("ecef", 3)

        for i in (1, 2, 3):
            _write_var(f, f"x{i}", (f"x{i}ghost",), xg[f"x{i}"])
            _write_var(f, f"x{i}i", (f"x{i}i",), xg[f"x{i}i"])
            _write_var(f, f"dx{i}b", (f"x{i}d",), xg[f"dx{i}b"])
            _write_var(f, f"dx{i}h", (f"x{i}",), xg[f"dx{i}h"])
            _write_var(f, f"h{i}", ("x3ghost", "x2ghost", "x1ghost"), xg[f"h{i}"].transpose())
            _write_var(f, f"h{i}x1i", ("x3", "x2", "x1i"), xg[f"h{i}x1i"].transpose())
            _write_var(f, f"h{i}x2i", ("x3", "x2i", "x1"), xg[f"h{i}x2i"].transpose())
            _write_var(f, f"h{i}x3i", ("x3i", "x2", "x1"), xg[f"h{i}x3i"].transpose())
            _write_var(f, f"gx{i}", ("x3", "x2", "x1"), xg[f"gx{i}"].transpose())
            _write_var(f, f"e{i}", ("ecef", "x3", "x2", "x1"), xg[f"e{i}"].transpose())

        for k in ("alt", "glat", "glon", "Bmag", "nullpts", "r", "theta", "phi", "x", "y", "z"):
            _write_var(f, k, ("x3", "x2", "x1"), xg[k].transpose())

        for k in ("er", "etheta", "ephi"):
            _write_var(f, k, ("ecef", "x3", "x2", "x1"), xg[k].transpose())

        _write_var(f, "I", ("x3", "x2"), xg["I"].transpose())


def _write_var(f, name: str, dims: tuple, value: np.ndarray):
    g = f.createVariable(name, np.float32, dims, zlib=True, complevel=1, shuffle=True, fletcher32=True, fill_value=np.nan)
    g[:] = value


def write_state(time: datetime, ns: np.ndarray, vs: np.ndarray, Ts: np.ndarray, out_dir: Path):
    """
     WRITE STATE VARIABLE DATA.
    NOTE THAT WE don't write ANY OF THE ELECTRODYNAMIC
    VARIABLES SINCE THEY ARE NOT NEEDED TO START THINGS
    UP IN THE FORTRAN CODE.

    INPUT ARRAYS SHOULD BE TRIMMED TO THE CORRECT SIZE
    I.E. THEY SHOULD NOT INCLUDE GHOST CELLS
    """

    fn = out_dir / "inputs/initial_conditions.nc"
    print("write_state:", fn)

    with Dataset(fn, "w") as f:
        p4 = (0, 3, 2, 1)

        f.createDimension("ymd", 3)
        g = f.createVariable("ymd", np.int32, "ymd")
        g[:] = [time.year, time.month, time.day]

        g = f.createVariable("UTsec", np.float32)
        g[:] = time.hour * 3600 + time.minute * 60 + time.second + time.microsecond / 1e6

        f.createDimension("species", 7)
        f.createDimension("x1", ns.shape[1])
        f.createDimension("x2", ns.shape[2])
        f.createDimension("x3", ns.shape[3])

        _write_var(f, "ns", ("species", "x3", "x2", "x1"), ns.transpose(p4))
        _write_var(f, "vsx1", ("species", "x3", "x2", "x1"), vs.transpose(p4))
        _write_var(f, "Ts", ("species", "x3", "x2", "x1"), Ts.transpose(p4))


def write_Efield(p: T.Dict[str, T.Any], E: T.Dict[str, np.ndarray]):
    """
    write Efield to disk

    TODO: verify dimensions vs. data vs. Fortran order
    """

    outdir = E["Efield_outdir"]
    print("write E-field data to", outdir)

    with Dataset(outdir / "simsize.nc", "w") as f:
        for k in ("llon", "llat"):
            g = f.createVariable(k, np.int32)
            g[:] = E[k]

    with Dataset(outdir / "simgrid.nc", "w") as f:
        f.createDimension("lon", E["mlon"].size)
        f.createDimension("lat", E["mlat"].size)
        _write_var(f, "mlon", ("lon",), E["mlon"])
        _write_var(f, "mlat", ("lat",), E["mlat"])

    for i, t in enumerate(E["time"]):
        fn = outdir / (datetime2ymd_hourdec(t) + ".nc")

        # FOR EACH FRAME WRITE A BC TYPE AND THEN OUTPUT BACKGROUND AND BCs
        with Dataset(fn, "w") as f:
            f.createDimension("lon", E["mlon"].size)
            f.createDimension("lat", E["mlat"].size)

            g = f.createVariable("flagdirich", np.int32)
            g[:] = p["flagdirich"]

            f.createDimension("ymd", 3)
            g = f.createVariable("ymd", np.int32, "ymd")
            g[:] = [t.year, t.month, t.day]

            g = f.createVariable("UTsec", np.float32)
            g[:] = t.hour * 3600 + t.minute * 60 + t.second + t.microsecond / 1e6

            for k in ("Exit", "Eyit", "Vminx1it", "Vmaxx1it"):
                _write_var(f, k, ("lat", "lon"), E[k][i, :, :].transpose())

            for k in ("Vminx2ist", "Vmaxx2ist"):
                _write_var(f, k, ("lat",), E[k][i, :])

            for k in ("Vminx3ist", "Vmaxx3ist"):
                _write_var(f, k, ("lon",), E[k][i, :])


def write_precip(precip: T.Dict[str, np.ndarray]):
    """
    write precipitation to disk

    TODO: verify dimensions vs. data vs. Fortran order
    """

    outdir = precip["precip_outdir"]
    print("write precipitation data to", outdir)

    with Dataset(outdir / "simsize.nc", "w") as f:
        for k in ("llon", "llat"):
            g = f.createVariable(k, np.int32)
            g[:] = precip[k]

    with Dataset(outdir / "simgrid.nc", "w") as f:
        f.createDimension("lon", precip["mlon"].size)
        f.createDimension("lat", precip["mlat"].size)
        _write_var(f, "mlon", ("lon",), precip["mlon"])
        _write_var(f, "mlat", ("lat",), precip["mlat"])

    for i, t in enumerate(precip["time"]):
        fn = outdir / (datetime2ymd_hourdec(t) + ".nc")

        with Dataset(fn, "w") as f:
            f.createDimension("lon", precip["mlon"].size)
            f.createDimension("lat", precip["mlat"].size)

            for k in ("Q", "E0"):
                _write_var(f, f"/{k}p", ("lat", "lon"), precip[k][i, :, :].transpose())
