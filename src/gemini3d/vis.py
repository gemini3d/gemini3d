from datetime import datetime
import numpy as np
import math
import argparse
from pathlib import Path
import scipy.interpolate as interp
import matplotlib as mpl
from matplotlib.pyplot import figure
from matplotlib.ticker import MultipleLocator
import typing as T

from .utils import gitrev
from .grid import readgrid
from .readdata import readdata

if T.TYPE_CHECKING:
    import matplotlib.axes as mpla
    import matplotlib.figure as mplf
    import matplotlib.collections as mplc

mpl.rcParams["axes.formatter.limits"] = (-3, 4)
mpl.rcParams["axes.formatter.useoffset"] = False
mpl.rcParams["axes.formatter.min_exponent"] = 4

R_EARTH = 6370e3
REF_ALT = 300  # km

CB_LBL = {
    "ne": "$n_e (m^{-3})$",
    "v1": "$v_1 (ms^{-1})$",
    "Ti": "$T_i$ (K)",
    "Te": "$T_e$ (K)",
    "J1": "$J_1 (Am^{-2})$",
    "v2": "$v_2 (ms^{-1})$",
    "v3": "$v_3 (ms^{-1})$",
    "J2": "$J_2 (Am^{-2})$",
    "J3": "$J_3 (Am^{-2})$",
    "Phitop": r"$\Phi_{top}$ (V)",
    "Vmaxx1it": r"(V)",
}

PARAMS = ["ne", "v1", "Ti", "Te", "J1", "v2", "v3", "J2", "J3", "Phitop"]


def gemini_plot():
    p = argparse.ArgumentParser()
    p.add_argument("direc", help="directory to plot")
    p.add_argument("--mayavi", help="do 3D Mayavi plots", action="store_true")
    p.add_argument("-s", "--saveplots", help="save plots to data directory", action="store_true")
    p.add_argument("--only", help="only plot these quantities", nargs="+")
    p = p.parse_args()

    if p.mayavi:
        from . import vis3d
        from mayavi.mlab import show
    else:
        from matplotlib.pyplot import show

    direc = Path(p.direc).expanduser().resolve(strict=True)
    if p.saveplots:
        from matplotlib.figure import Figure

        fg = Figure(constrained_layout=True)
        save_dir = direc / "plots"
        save_dir.mkdir(parents=True, exist_ok=True)
    else:
        fg = None
        save_dir = None

    grid = readgrid(direc)

    flist = sorted(direc.glob("*.h5"))
    if len(flist) == 0:
        flist = sorted(direc.glob("*.dat"))
    # %% loop over files / time
    for file in flist:
        try:
            dat = readdata(file)
        except Exception as e:
            print(f"SKIP: {file}   {e}")
            continue
        if "mlon" in dat and "mlon" not in grid:
            grid["mlon"] = dat["mlon"]
            grid["mlat"] = dat["mlat"]

        if p.mayavi:
            vis3d.plotframe(grid, dat, params=p.only, save_dir=save_dir)
        else:
            plotframe(grid, dat, params=p.only, save_dir=save_dir, fg=fg)

        if p.saveplots:
            print(f"{dat['time']} => {save_dir}")
        else:
            show()


def plotframe(
    grid: T.Dict[str, np.ndarray],
    dat: T.Dict[str, T.Any],
    params: T.Sequence[str] = None,
    save_dir: Path = None,
    fg: "mplf.Figure" = None,
):
    """
    if save_dir, plots will not be visible while generating to speed plot writing
    """
    if not params:
        params = PARAMS

    plotfun = grid2plotfun(grid)

    time = dat["time"]

    for k in params:
        if k not in dat:  # not present at this time step, often just the first time step
            continue
        if save_dir is None or fg is None:
            fg = figure(num=k, constrained_layout=True)
        fg.clf()
        plotfun(time, grid, dat[k][1].squeeze(), k, fg, wavelength=dat.get("wavelength"))
        if save_dir:
            fg.savefig(save_dir / f"{k}-{time.isoformat().replace(':','')}.png")


def grid2plotfun(grid: T.Dict[str, np.ndarray]):
    plotfun = None
    h1 = grid.get("h1")
    for k in ("lx", "lxs", "lx1"):
        if k in grid:
            if k == "lx1":
                lxs = (grid["lx1"], grid["lx2"], grid["lx3"])
            else:
                lxs = grid[k]
    if h1 is not None:
        minh1 = h1.min()
        maxh1 = h1.max()
        if (abs(minh1 - 1) > 1e-4) or (abs(maxh1 - 1) > 1e-4):  # curvilinear grid
            if (lxs[1] > 1) and (lxs[2] > 1):
                plotfun = plot3D_curv_frames_long
            else:
                plotfun = plot2D_curv
    if plotfun is None:  # cartesian grid
        if (lxs[1] > 1) and (lxs[2] > 1):
            plotfun = plot3D_cart_frames_long_ENU
        else:
            plotfun = plot2D_cart

    return plotfun


def plot3D_curv_frames_long(
    time: datetime,
    grid: T.Dict[str, np.ndarray],
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
    **kwargs,
):
    raise NotImplementedError


def plot2D_curv(
    time: datetime,
    grid: T.Dict[str, np.ndarray],
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
    **kwargs,
):
    raise NotImplementedError


def plot12(
    x: np.ndarray,
    z: np.ndarray,
    parm: np.ndarray,
    name: str,
    cmap: str,
    vmin: float,
    vmax: float,
    fg: "mplf.Figure",
    ax: "mpla.Axes" = None,
) -> "mplc.QuadMesh":

    if parm.ndim != 2:
        raise ValueError(f"data must have 2 dimensions, you have {parm.shape}")

    if ax is None:
        ax = fg.gca()

    hi = ax.pcolormesh(x / 1e3, z / 1e3, parm, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.set_xlabel("eastward dist. (km)")
    ax.set_ylabel("upward dist. (km)")
    ax.axhline(REF_ALT, color="w", linestyle="--", linewidth=2)
    fg.colorbar(hi, ax=ax, label=CB_LBL[name])

    return hi


def plot13(
    y: np.ndarray,
    z: np.ndarray,
    parm: np.ndarray,
    name: str,
    cmap: str,
    vmin: float,
    vmax: float,
    fg: "mplf.Figure",
    ax: "mpla.Axes" = None,
) -> "mplc.QuadMesh":

    if parm.ndim != 2:
        raise ValueError(f"data must have 2 dimensions, you have {parm.shape}")

    if ax is None:
        ax = fg.gca()

    hi = ax.pcolormesh(y / 1e3, z / 1e3, parm, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.set_xlabel("northward dist. (km)")
    ax.set_ylabel("upward dist. (km)")
    fg.colorbar(hi, ax=ax, label=CB_LBL[name])

    return hi


def plot23(
    x: np.ndarray,
    y: np.ndarray,
    parm: np.ndarray,
    name: str,
    cmap: str,
    vmin: float,
    vmax: float,
    fg: "mplf.Figure",
    ax: "mpla.Axes" = None,
) -> "mplc.QuadMesh":

    if parm.ndim != 2:
        raise ValueError(f"data must have 2 dimensions, you have {parm.shape}")

    if ax is None:
        ax = fg.gca()

    hi = ax.pcolormesh(x / 1e3, y / 1e3, parm, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlabel("eastward dist. (km)")
    ax.set_ylabel("northward dist. (km)")
    fg.colorbar(hi, ax=ax, label=CB_LBL[name])

    return hi


def plot1d2(x: np.ndarray, parm: np.ndarray, name: str, fg: "mplf.Figure", ax: "mpla.Axes" = None):

    if parm.ndim != 1:
        raise ValueError("expecting 1-D data oriented east-west (along latitude)")

    if ax is None:
        ax = fg.gca()

    ax.plot(x / 1e3, parm)
    ax.set_xlabel("eastward dist. (km)")
    ax.set_ylabel(CB_LBL[name])


def plot1d3(y: np.ndarray, parm: np.ndarray, name: str, fg: "mplf.Figure", ax: "mpla.Axes" = None):

    if parm.ndim != 1:
        raise ValueError("expecting 1-D data oriented east-west (along latitude)")

    if ax is None:
        ax = fg.gca()

    ax.plot(y / 1e3, parm)
    ax.set_xlabel("northward dist. (km)")
    ax.set_ylabel(CB_LBL[name])


def plot_interp(
    time: datetime,
    grid: T.Dict[str, np.ndarray],
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
    **kwargs,
):
    """

    xp:  eastward distance (rads.)
        should be interpreted as northward distance (in rads.).
        Irrespective of ordering of xg.theta, this will be monotonic increasing!!!

    zp: altitude (meters)

    y:  this is a mag colat. coordinate and is only used for defining
        grid in linspaces below, runs backward from north distance,
        hence the negative sign
    """

    cmap = None
    is_Efield = False
    vmin = None
    vmax = None
    if name.startswith("J") or name == "Phitop":
        cmap = "bwr"
        vmax = abs(parm).max()
        vmin = -vmax
    elif name.startswith("v"):
        cmap = "bwr"
        vmax = 80.0
        vmin = -vmax
    elif name.startswith(("V", "E")):
        is_Efield = True
        cmap = "bwr"
        vmax = abs(parm).max()
        vmin = -vmax
    elif name.startswith("T"):
        vmin = 0.0
        vmax = parm.max()
    elif name.startswith("n"):
        vmin = 1e-7

    # %% SIZE OF SIMULATION
    for k in ("lx", "lxs", "lx1"):
        if k in grid:
            if k == "lx1":
                lxs = (grid["lx1"], grid["lx2"], grid["lx3"])
            else:
                lxs = grid[k]
    lx1, lx2, lx3 = lxs
    inds1 = slice(2, lx1 + 2)
    inds2 = slice(2, lx2 + 2)
    inds3 = slice(2, lx3 + 2)
    # %% SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
    meantheta = grid["theta"].mean()
    # this is a mag colat. coordinate and is only used for defining grid in linspaces below
    # runs backward from north distance, hence the negative sign
    # [radians]
    y = -1 * (grid["theta"] - meantheta)
    # eastward distance [radians]
    x = grid["x2"][inds2] / R_EARTH / math.sin(meantheta)
    # altitude [meters]
    z = grid["alt"] / 1e3

    # arbitrary output plot resolution
    lxp = 500
    lyp = 500
    lzp = 500

    # eastward distance [meters]
    xp = np.linspace(x.min(), x.max(), lxp) * R_EARTH * math.sin(meantheta)
    # northward distance [meters]
    yp = np.linspace(y.min(), y.max(), lyp) * R_EARTH
    # upward distance [meters]
    zp = np.linspace(z.min(), z.max(), lzp) * 1e3

    # %% INTERPOLATE ONTO PLOTTING GRID
    if lxs[2] == 1:  # alt./lon. slice
        ax = fg.gca()
        ax.set_title(f"{name}: {time.isoformat()}  {gitrev()}")
        # meridional meshgrid, this defines the grid for plotting
        # slice expects the first dim. to be "y" ("z" in the 2D case)
        # %% CONVERT ANGULAR COORDINATES TO MLAT,MLON
        i = np.argsort(xp)  # FIXME: this was in Matlab code--what is its purpose?

        if name == "rayleighs":
            f = interp.interp1d(grid["x2"][inds2], parm, axis=1, bounds_error=False)
            # hack for pcolormesh to put labels in center of pixel
            wl = kwargs["wavelength"] + [""]
            hi = ax.pcolormesh(xp / 1e3, np.arange(len(wl)), f(xp)[:, i])
            ax.set_yticks(np.arange(len(wl)) + 0.5)
            ax.set_yticklabels(wl)
            ax.set_ylim(0, len(wl) - 1)
            # end hack
            ax.set_ylabel(r"wavelength $\AA$")
            ax.set_xlabel("eastward dist. (km)")
        elif parm.ndim == 2:
            f = interp.interp2d(grid["x2"][inds2], grid["x1"][inds1], parm, bounds_error=False)
            plot12(xp[i], zp, f(xp, zp)[:, i], name, cmap, vmin, vmax, fg, ax)
        elif parm.ndim == 1:  # phitop
            f = interp.interp1d(grid["x2"][inds2], parm, bounds_error=False)
            plot1d2(xp, f(xp), name, fg, ax)
        else:
            raise ValueError(f"{name}: only 2D and 1D data are expected--squeeze data")
    elif lxs[1] == 1:  # alt./lat. slice
        ax = fg.gca()
        ax.set_title(f"{name}: {time.isoformat()}  {gitrev()}")
        # so north dist, east dist., alt.
        # slice expects the first dim. to be "y"
        # %% CONVERT ANGULAR COORDINATES TO MLAT,MLON
        i = np.argsort(yp)  # FIXME: this was in Matlab code--what is its purpose?

        if name == "rayleighs":
            # FIXME: this needs to be tested
            f = interp.interp1d(grid["x3"][inds3], parm, axis=1, bounds_error=False)
            # hack for pcolormesh to put labels in center of pixel
            wl = kwargs["wavelength"] + [""]
            hi = ax.pcolormesh(np.arange(len(wl)), yp / 1e3, f(yp)[:, i].T)
            ax.set_xticks(np.arange(len(wl)) + 0.5)
            ax.set_xticklabels(wl)
            ax.set_xlim(0, len(wl) - 1)
            # end hack
            ax.set_xlabel(r"wavelength $\AA$")
            ax.set_ylabel("northward dist. (km)")
        elif parm.ndim == 2:
            f = interp.interp2d(grid["x3"][inds3], grid["x1"][inds1], parm, bounds_error=False)
            parmp = f(yp, zp).reshape((lzp, lyp))
            plot13(yp[i], zp, parmp[:, i], name, cmap, vmin, vmax, fg, ax)
        elif parm.ndim == 1:  # phitop
            f = interp.interp1d(grid["x3"][inds3], parm, bounds_error=False)
            plot1d3(yp, f(yp), name, fg, ax)
        else:
            raise ValueError(f"{name}: only 2D and 1D data are expected--squeeze data")

    else:  # 3-panel plot, vs. single-panel plots of 2-D cases
        if name == "rayleighs":
            axs = fg.subplots(2, 2, sharey=True, sharex=True).ravel()
            fg.suptitle(f"{name}: {time.isoformat()}  {gitrev()}", y=0.99)
            # arbitrary pick of which emission lines to plot lat/lon slices
            for j, i in enumerate([1, 3, 4, 8]):
                f = interp.interp2d(
                    grid["x3"][inds3], grid["x2"][inds2], parm[i, :, :], bounds_error=False
                )
                hi = axs[j].pcolormesh(xp / 1e3, yp / 1e3, f(yp, xp))
                axs[j].set_title(kwargs["wavelength"][i] + r"$\AA$")
                fg.colorbar(hi, ax=axs[j], label="Rayleighs")
            axs[2].set_xlabel("eastward dist. (km)")
            axs[2].set_ylabel("northward dist. (km)")
            return
        elif parm.ndim == 3:
            fg.set_size_inches((18, 5))
            axs = fg.subplots(1, 3, sharey=False, sharex=False)
            fg.suptitle(f"{name}: {time.isoformat()}  {gitrev()}", y=0.99)
        elif is_Efield:
            # like phitop, SINGLE plot
            ax = fg.gca()
            hi = ax.pcolormesh(grid["mlon"], grid["mlat"], parm, cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_xlabel("magnetic longitude (deg.)")
            ax.set_ylabel("magnetic latitude (deg.)")
            ax.set_title(f"{name}: {time.isoformat()}  {gitrev()}")
            fg.colorbar(hi, ax=ax, label=CB_LBL[name])
            return
        else:
            # like phitop, SINGLE plot
            ax = fg.gca()
            f = interp.interp2d(grid["x3"][inds3], grid["x2"][inds2], parm, bounds_error=False)
            hi = ax.pcolormesh(xp / 1e3, yp / 1e3, f(yp, xp), cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_xlabel("eastward dist. (km)")
            ax.set_ylabel("northward dist. (km)")
            ax.set_title(f"{name}: {time.isoformat()}  {gitrev()}")
            fg.colorbar(hi, ax=ax, label=CB_LBL[name])
            return

        # %% CONVERT TO DISTANCE UP, EAST, NORTH (left panel)
        # JUST PICK AN X3 LOCATION FOR THE MERIDIONAL SLICE PLOT,
        # AND AN ALTITUDE FOR THE LAT./LON. SLICE
        ix3 = lx3 // 2 - 1  # arbitrary slice, to match Matlab
        f = interp.interp2d(
            grid["x2"][inds2], grid["x1"][inds1], parm[:, :, ix3], bounds_error=False
        )
        # CONVERT ANGULAR COORDINATES TO MLAT,MLON
        ix = np.argsort(xp)
        iy = np.argsort(yp)
        plot12(xp[ix], zp, f(xp, zp)[:, ix], name, cmap, vmin, vmax, fg, axs[0])
        # %% LAT./LONG. SLICE COORDINATES (center panel)
        zp2 = REF_ALT
        X3, Y3, Z3 = np.meshgrid(xp, yp, zp2 * 1e3)
        # transpose: so north dist, east dist., alt.
        parmp = interp.interpn(
            points=(grid["x1"][inds1], grid["x2"][inds2], grid["x3"][inds3]),
            values=parm,
            xi=np.column_stack((Z3.ravel(), X3.ravel(), Y3.ravel())),
            bounds_error=False,
        ).reshape((1, lxp, lyp))

        parmp = parmp[:, :, iy]  # must be indexed in two steps
        plot23(xp[ix], yp[iy], parmp[0, ix, :], name, cmap, vmin, vmax, fg, axs[1])
        # %% ALT/LAT SLICE (right panel)
        ix2 = lx2 // 2 - 1  # arbitrary slice, to match Matlab
        f = interp.interp2d(
            grid["x3"][inds3], grid["x1"][inds1], parm[:, ix2, :], bounds_error=False
        )
        hi = plot13(yp[iy], zp, f(yp, zp)[:, iy], name, cmap, vmin, vmax, fg, axs[2])

        fg.colorbar(hi, ax=axs, aspect=60, pad=0.01)


plot3D_cart_frames_long_ENU = plot_interp
plot2D_cart = plot_interp
