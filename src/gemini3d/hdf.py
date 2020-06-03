from pathlib import Path
import typing as T
import numpy as np
import logging
import h5py
from datetime import datetime

from .utils import datetime2ymd_hourdec, ymdhourdec2datetime

LSP = 7


def get_simsize(path: Path) -> T.Tuple[int, ...]:
    """
    get simulation size
    """
    path = Path(path).expanduser().resolve()

    with h5py.File(path, "r") as f:
        if "lxs" in f:
            lxs = f["lxs"][:]
        elif "lx" in f:
            lxs = f["lx"][:]
        elif "lx1" in f:
            if f["lx1"].ndim > 0:
                lxs = (
                    f["lx1"][:].squeeze()[()],
                    f["lx2"][:].squeeze()[()],
                    f["lx3"][:].squeeze()[()],
                )
            else:
                lxs = (f["lx1"][()], f["lx2"][()], f["lx3"][()])
        else:
            raise KeyError(f"could not find '/lxs', '/lx' or '/lx1' in {path.as_posix()}")

    return lxs


def read_state(fn: Path) -> T.Dict[str, T.Any]:
    """
    load initial condition data
    """

    with h5py.File(fn, "r") as f:
        return {"ns": f["/ns"][:], "vs": f["/vsx1"][:], "Ts": f["/Ts"][:]}


def write_state(time: datetime, ns: np.ndarray, vs: np.ndarray, Ts: np.ndarray, fn: Path):
    """
    write STATE VARIABLE initial conditions

    NOTE THAT WE don't write ANY OF THE ELECTRODYNAMIC
    VARIABLES SINCE THEY ARE NOT NEEDED TO START THINGS
    UP IN THE FORTRAN CODE.

    INPUT ARRAYS SHOULD BE TRIMMED TO THE CORRECT SIZE
    I.E. THEY SHOULD NOT INCLUDE GHOST CELLS

    NOTE: The .transpose() reverses the dimension order.
    The HDF Group never implemented the intended H5T_array_create(..., perm)
    and it's deprecated.
    Fortran, including the HDF Group Fortran interfaces and h5fortran as well as
    Matlab read/write HDF5 in Fortran order. h5py read/write HDF5 in C order so we
    need the .transpose() for h5py
    """

    print("hdf:write_state:", fn)

    with h5py.File(fn, "w") as f:
        f["/ymd"] = [time.year, time.month, time.day]
        f["/UTsec"] = time.hour * 3600 + time.minute * 60 + time.second + time.microsecond / 1e6

        p4 = (0, 3, 2, 1)
        # we have to reverse axes order and put lsp at the last dim

        f.create_dataset(
            "/ns",
            data=ns.transpose(p4),
            dtype=np.float32,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
            fletcher32=True,
        )
        f.create_dataset(
            "/vsx1",
            data=vs.transpose(p4),
            dtype=np.float32,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
            fletcher32=True,
        )
        f.create_dataset(
            "/Ts",
            data=Ts.transpose(p4),
            dtype=np.float32,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
            fletcher32=True,
        )


def readgrid(fn: Path) -> T.Dict[str, np.ndarray]:
    """
    get simulation grid

    Parameters
    ----------
    fn: pathlib.Path
        filepath to simgrid.h5

    Returns
    -------
    grid: dict
        grid parameters
    """

    grid: T.Dict[str, T.Any] = {}

    if not fn.is_file():
        logging.error(f"{fn} grid file is not present. Will try to load rest of data.")
        return grid

    grid["lxs"] = get_simsize(fn.with_name("simsize.h5"))

    with h5py.File(fn, "r") as f:
        for key in f.keys():
            grid[key] = f[key][:]

    return grid


def write_grid(size_fn: Path, grid_fn: Path, xg: T.Dict[str, T.Any]):
    """ writes grid to disk

    Parameters
    ----------

    size_fn: pathlib.Path
        file to write
    grid_fn: pathlib.Path
        file to write
    xg: dict
        grid values

    NOTE: The .transpose() reverses the dimension order.
    The HDF Group never implemented the intended H5T_array_create(..., perm)
    and it's deprecated.
    Fortran, including the HDF Group Fortran interfaces and h5fortran as well as
    Matlab read/write HDF5 in Fortran order. h5py read/write HDF5 in C order so we
    need the .transpose() for h5py
    """

    print("hdf:write_grid:", size_fn)
    with h5py.File(size_fn, "w") as h:
        h["/lx"] = xg["lx"]

    print("hdf:write_grid:", grid_fn)
    with h5py.File(grid_fn, "w") as h:
        for i in (1, 2, 3):
            for k in (
                f"x{i}",
                f"x{i}i",
                f"dx{i}b",
                f"dx{i}h",
                f"h{i}",
                f"h{i}x1i",
                f"h{i}x2i",
                f"h{i}x3i",
                f"gx{i}",
                f"e{i}",
            ):
                if xg[k].ndim >= 2:
                    h.create_dataset(
                        f"/{k}",
                        data=xg[k].transpose(),
                        dtype=np.float32,
                        compression="gzip",
                        compression_opts=1,
                        shuffle=True,
                        fletcher32=True,
                    )
                else:
                    h[f"/{k}"] = xg[k].astype(np.float32)

        for k in (
            "alt",
            "glat",
            "glon",
            "Bmag",
            "I",
            "nullpts",
            "er",
            "etheta",
            "ephi",
            "r",
            "theta",
            "phi",
            "x",
            "y",
            "z",
        ):
            if xg[k].ndim >= 2:
                h.create_dataset(
                    f"/{k}",
                    data=xg[k].transpose(),
                    dtype=np.float32,
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                    fletcher32=True,
                )
            else:
                h[f"/{k}"] = xg[k].astype(np.float32)


def read_Efield(fn: Path) -> T.Dict[str, T.Any]:
    """
    load electric field
    """

    # sizefn = fn.with_name("simsize.h5")  # NOT the whole sim simsize.dat
    # with h5py.File(sizefn, "r") as f:
    #     E["llon"] = f["/llon"][()]
    #     E["llat"] = f["/llat"][()]

    gridfn = fn.with_name("simgrid.h5")  # NOT the whole sim simgrid.dat
    with h5py.File(gridfn, "r") as f:
        E = {"mlon": f["/mlon"][:], "mlat": f["/mlat"][:]}

    with h5py.File(fn, "r") as f:
        E["flagdirich"] = f["flagdirich"]
        for p in ("Exit", "Eyit", "Vminx1it", "Vmaxx1it"):
            E[p] = [("x2", "x3"), f[p][:]]
        for p in ("Vminx2ist", "Vmaxx2ist"):
            E[p] = [("x2",), f[p][:]]
        for p in ("Vminx3ist", "Vmaxx3ist"):
            E[p] = [("x3",), f[p][:]]

    return E


def write_Efield(outdir: Path, E: T.Dict[str, np.ndarray]):
    """
    write Efield to disk
    """

    with h5py.File(outdir / "simsize.h5", "w") as f:
        f["/llon"] = E["llon"]
        f["/llat"] = E["llat"]

    with h5py.File(outdir / "simgrid.h5", "w") as f:
        f["/mlon"] = E["mlon"].astype(np.float32)
        f["/mlat"] = E["mlat"].astype(np.float32)

    for i, t in enumerate(E["time"]):
        fn = outdir / (datetime2ymd_hourdec(t) + ".h5")

        # FOR EACH FRAME WRITE A BC TYPE AND THEN OUTPUT BACKGROUND AND BCs
        with h5py.File(fn, "w") as f:
            f["/flagdirich"] = E["flagdirich"][i].astype(np.int32)
            f["/time/ymd"] = [t.year, t.month, t.day]
            f["/time/UTsec"] = t.hour * 3600 + t.minute * 60 + t.second + t.microsecond / 1e6

            for k in ("Exit", "Eyit", "Vminx1it", "Vmaxx1it"):
                f.create_dataset(
                    f"/{k}",
                    data=E[k][i, :, :].transpose(),
                    dtype=np.float32,
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                    fletcher32=True,
                )
            for k in ("Vminx2ist", "Vmaxx2ist", "Vminx3ist", "Vmaxx3ist"):
                f[f"/{k}"] = E[k][i, :].astype(np.float32)


def read_precip(fn: Path) -> T.Dict[str, T.Any]:

    # with h5py.File(path / "simsize.h5", "r") as f:
    #     dat["llon"] = f["/llon"][()]
    #     dat["llat"] = f["/llat"][()]

    with h5py.File(fn.with_name("simgrid.h5"), "r") as f:
        dat = {"mlon": f["/mlon"][:], "mlat": f["/mlat"][:]}

    with h5py.File(fn, "r") as f:
        for k in ("Q", "E0"):
            dat[k] = f[f"/{k}p"][:]

    return dat


def write_precip(outdir: Path, precip: T.Dict[str, T.Any]):

    with h5py.File(outdir / "simsize.h5", "w") as f:
        f["/llon"] = precip["llon"]
        f["/llat"] = precip["llat"]

    with h5py.File(outdir / "simgrid.h5", "w") as f:
        f["/mlon"] = precip["mlon"].astype(np.float32)
        f["/mlat"] = precip["mlat"].astype(np.float32)

    for i, t in enumerate(precip["time"]):
        fn = outdir / (datetime2ymd_hourdec(t) + ".h5")

        with h5py.File(fn, "w") as f:
            for k in ("Q", "E0"):
                f.create_dataset(
                    f"/{k}p",
                    data=precip[k][i, :, :].transpose(),
                    dtype=np.float32,
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                    fletcher32=True,
                )


def loadframe3d_curv(fn: Path, lxs: T.Sequence[int]) -> T.Dict[str, T.Any]:
    """
    end users should normally use loadframe() instead
    """

    #    grid = readgrid(fn.parent / "inputs/simgrid.h5")
    #    dat = xarray.Dataset(
    #        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    #    )

    dat: T.Dict[str, T.Any] = {}

    with h5py.File(fn, "r") as f:
        dat["time"] = ymdhourdec2datetime(
            f["time/ymd"][0], f["time/ymd"][1], f["time/ymd"][2], f["time/UThour"][()]
        )

        if lxs[2] == 1:  # east-west
            p4 = (0, 3, 1, 2)
            p3 = (2, 0, 1)
        else:  # 3D or north-south, no swap
            p4 = (0, 3, 2, 1)
            p3 = (2, 1, 0)

        ns = f["/nsall"][:].transpose(p4)
        # np.any() in case neither is an np.ndarray
        if ns.shape[0] != 7 or np.any(ns.shape[1:] != lxs):
            raise ValueError(
                f"may have wrong permutation on read. lxs: {lxs}  ns x1,x2,x3: {ns.shape}"
            )
        dat["ns"] = (("lsp", "x1", "x2", "x3"), ns)
        vs = f["/vs1all"][:].transpose(p4)
        dat["vs"] = (("lsp", "x1", "x2", "x3"), vs)
        Ts = f["/Tsall"][:].transpose(p4)
        dat["Ts"] = (("lsp", "x1", "x2", "x3"), Ts)

        dat["ne"] = (("x1", "x2", "x3"), ns[LSP - 1, :, :, :])

        dat["v1"] = (
            ("x1", "x2", "x3"),
            (ns[:6, :, :, :] * vs[:6, :, :, :]).sum(axis=0) / dat["ne"][1],
        )

        dat["Ti"] = (
            ("x1", "x2", "x3"),
            (ns[:6, :, :, :] * Ts[:6, :, :, :]).sum(axis=0) / dat["ne"][1],
        )
        dat["Te"] = (("x1", "x2", "x3"), Ts[LSP - 1, :, :, :])

        dat["J1"] = (("x1", "x2", "x3"), f["/J1all"][:].transpose(p3))
        # np.any() in case neither is an np.ndarray
        if np.any(dat["J1"][1].shape != lxs):
            raise ValueError("may have wrong permutation on read")
        dat["J2"] = (("x1", "x2", "x3"), f["/J2all"][:].transpose(p3))
        dat["J3"] = (("x1", "x2", "x3"), f["/J3all"][:].transpose(p3))

        dat["v2"] = (("x1", "x2", "x3"), f["/v2avgall"][:].transpose(p3))
        dat["v3"] = (("x1", "x2", "x3"), f["/v3avgall"][:].transpose(p3))

        dat["Phitop"] = (("x2", "x3"), f["/Phiall"][:].transpose())

    return dat


def loadframe3d_curvavg(fn: Path, lxs: T.Sequence[int]) -> T.Dict[str, T.Any]:
    """
    end users should normally use loadframe() instead

    Parameters
    ----------
    fn: pathlib.Path
        filename of this timestep of simulation output
    """
    #    grid = readgrid(fn.parent / "inputs/simgrid.h5")
    #    dat = xarray.Dataset(
    #        coords={"x1": grid["x1"][2:-2], "x2": grid["x2"][2:-2], "x3": grid["x3"][2:-2]}
    #    )
    dat: T.Dict[str, T.Any] = {}

    with h5py.File(fn, "r") as f:
        dat["time"] = ymdhourdec2datetime(
            f["time/ymd"][0], f["time/ymd"][1], f["time/ymd"][2], f["/time/UThour"][()]
        )

        dat["ne"] = [("x1", "x2", "x3"), f["/neall"][:].transpose(2, 0, 1)]
        dat["v1"] = [("x1", "x2", "x3"), f["/v1avgall"][:].transpose(2, 0, 1)]
        dat["Ti"] = [("x1", "x2", "x3"), f["/Tavgall"][:].transpose(2, 0, 1)]
        dat["Te"] = [("x1", "x2", "x3"), f["/TEall"][:].transpose(2, 0, 1)]
        dat["J1"] = [("x1", "x2", "x3"), f["/J1all"][:].transpose(2, 0, 1)]
        dat["J2"] = [("x1", "x2", "x3"), f["/J2all"][:].transpose(2, 0, 1)]
        dat["J3"] = [("x1", "x2", "x3"), f["/J3all"][:].transpose(2, 0, 1)]
        dat["v2"] = [("x1", "x2", "x3"), f["/v2avgall"][:].transpose(2, 0, 1)]
        dat["v3"] = [("x1", "x2", "x3"), f["/v3avgall"][:].transpose(2, 0, 1)]
        dat["Phitop"] = [("x2", "x3"), f["/Phiall"][:]]

    return dat


def loadglow_aurmap(fn: Path) -> T.Dict[str, T.Any]:
    """
    read the auroral output from GLOW

    Parameters
    ----------
    fn: pathlib.Path
        filename of this timestep of simulation output
    """

    with h5py.File(fn, "r") as h:
        dat = {"rayleighs": [("wavelength", "x2", "x3"), h["/aurora/iverout"][:]]}

    return dat
