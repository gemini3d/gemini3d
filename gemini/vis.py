from datetime import datetime
import numpy as np
import math
from pathlib import Path
import scipy.interpolate as interp
import matplotlib as mpl
from matplotlib.pyplot import figure
from matplotlib.ticker import MultipleLocator
import typing

from .utils import gitrev

if typing.TYPE_CHECKING:
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
}


def plotframe(
    grid: typing.Dict[str, np.ndarray], dat: typing.Dict[str, typing.Any], save_dir: Path = None, fg: "mplf.Figure" = None
):
    """
    if save_dir, plots will not be visible while generating to speed plot writing
    """
    plotfun = grid2plotfun(grid)

    time = dat["time"]

    for k in ("ne", "v1", "Ti", "Te", "J1", "v2", "v3", "J2", "J3", "Phitop"):
        if save_dir is None or fg is None:
            fg = figure(num=k, constrained_layout=True)
        fg.clf()
        plotfun(time, grid, dat[k][1].squeeze(), k, fg)
        if save_dir:
            fg.savefig(save_dir / f"{k}-{time.isoformat().replace(':','')}.png")


def grid2plotfun(grid: typing.Dict[str, np.ndarray]):
    plotfun = None
    h1 = grid.get("h1")
    if h1 is not None:
        minh1 = h1.min()
        maxh1 = h1.max()
        if (abs(minh1 - 1) > 1e-4) or (abs(maxh1 - 1) > 1e-4):  # curvilinear grid
            if (grid["lx"][1] > 1) and (grid["lx"][2] > 1):
                plotfun = plot3D_curv_frames_long
            else:
                plotfun = plot2D_curv
    if plotfun is None:  # cartesian grid
        if (grid["lx"][1] > 1) and (grid["lx"][2] > 1):
            plotfun = plot3D_cart_frames_long_ENU
        else:
            plotfun = plot2D_cart

    return plotfun


def plot3D_curv_frames_long(time: datetime, grid: typing.Dict[str, np.ndarray], parm: np.ndarray, name: str, fg: "mplf.Figure"):
    raise NotImplementedError


def plot2D_curv(time: datetime, grid: typing.Dict[str, np.ndarray], parm: np.ndarray, name: str, fg: "mplf.Figure"):
    raise NotImplementedError


def plot3D_cart_frames_long_ENU(time: datetime, grid: typing.Dict[str, np.ndarray], parm: np.ndarray, name: str, fg: "mplf.Figure"):

    plot_interp(time, grid, parm, name, fg)


def plot2D_cart(time: datetime, grid: typing.Dict[str, np.ndarray], parm: np.ndarray, name: str, fg: "mplf.Figure"):

    plot_interp(time, grid, parm, name, fg)


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
        make_colorbar = True
    else:
        make_colorbar = False

    hi = ax.pcolormesh(x / 1e3, z / 1e3, parm, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.set_xlabel("eastward dist. (km)")
    ax.set_ylabel("upward dist. (km)")
    ax.axhline(REF_ALT, color="w", linestyle="--", linewidth=2)
    if make_colorbar:
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
        make_colorbar = True
    else:
        make_colorbar = False

    hi = ax.pcolormesh(y / 1e3, z / 1e3, parm, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.set_xlabel("northward dist. (km)")
    ax.set_ylabel("upward dist. (km)")
    if make_colorbar:
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
        make_colorbar = True
    else:
        make_colorbar = False

    hi = ax.pcolormesh(x / 1e3, y / 1e3, parm, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlabel("eastward dist. (km)")
    ax.set_ylabel("northward dist. (km)")
    if make_colorbar:
        fg.colorbar(hi, ax=ax, label=CB_LBL[name])

    return hi


def plot1d2(x: np.ndarray, parm: np.ndarray, name: str, fg: "mplf.Figure", ax: "mpla.Axes" = None):

    if parm.ndim != 1:
        raise ValueError('expecting 1-D data oriented east-west (along latitude)')

    if ax is None:
        ax = fg.gca()

    ax.plot(x / 1e3, parm)
    ax.set_xlabel('eastward dist. (km)')
    ax.set_ylabel(CB_LBL[name])


def plot1d3(y: np.ndarray, parm: np.ndarray, name: str, fg: "mplf.Figure", ax: "mpla.Axes" = None):

    if parm.ndim != 1:
        raise ValueError('expecting 1-D data oriented east-west (along latitude)')

    if ax is None:
        ax = fg.gca()

    ax.plot(y / 1e3, parm)
    ax.set_xlabel('northward dist. (km)')
    ax.set_ylabel(CB_LBL[name])


def plot_interp(time: datetime, grid: typing.Dict[str, np.ndarray], parm: np.ndarray, name: str, fg: "mplf.Figure"):
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
    if name.startswith(("J", "v")) or name == 'Phitop':
        cmap = "bwr"
        vmax = abs(parm).max()
        vmin = -vmax
    elif name.startswith("T"):
        vmin = 0.0
        vmax = parm.max()
    else:
        vmin = vmax = None
    # %% SIZE OF SIMULATION
    lx1, lx2, lx3 = grid["lx"]
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
    if grid["lx"][2] == 1:  # alt./lon. slice
        ax = fg.gca()
        ax.set_title(f"{name}: {time.isoformat()}  {gitrev()}")
        # meridional meshgrid, this defines the grid for plotting
        # slice expects the first dim. to be "y" ("z" in the 2D case)
        # %% CONVERT ANGULAR COORDINATES TO MLAT,MLON
        i = np.argsort(xp)  # FIXME: this was in Matlab code--what is its purpose?

        if parm.ndim == 2:
            f = interp.interp2d(grid["x2"][inds2], grid["x1"][inds1], parm, bounds_error=False)
            parmp = f(xp, zp)
            plot12(xp[i], zp, parmp[:, i], name, cmap, vmin, vmax, fg)
        elif parm.ndim == 1:
            f = interp.interp1d(grid['x2'][inds2], parm, bounds_error=False)
            parmp = f(xp)
            plot1d2(xp, parmp, name, fg)
        else:
            raise ValueError(f'{name}: only 2D and 1D data are expected--squeeze data')
    elif grid["lx"][1] == 1:  # alt./lat. slice
        ax = fg.gca()
        ax.set_title(f"{name}: {time.isoformat()}  {gitrev()}")
        # so north dist, east dist., alt.
        # slice expects the first dim. to be "y"
        # %% CONVERT ANGULAR COORDINATES TO MLAT,MLON
        i = np.argsort(yp)  # FIXME: this was in Matlab code--what is its purpose?

        if parm.ndim == 2:
            f = interp.interp2d(grid["x3"][inds3], grid["x1"][inds1], parm, bounds_error=False)
            parmp = f(yp, zp).reshape((lzp, lyp))
            plot13(yp[i], zp, parmp[:, i], name, cmap, vmin, vmax, fg)
        elif parm.ndim == 1:
            f = interp.interp1d(grid['x3'][inds3], parm, bounds_error=False)
            parmp = f(yp)
            plot1d3(yp, parmp, name, fg)
        else:
            raise ValueError(f'{name}: only 2D and 1D data are expected--squeeze data')

    else:  # 3-panel plot, vs. single-panel plots of 2-D cases
        if parm.ndim == 3:
            fg.set_size_inches((18, 5))
            axs = fg.subplots(1, 3, sharey=False, sharex=False)
            fg.suptitle(f"{name}: {time.isoformat()}  {gitrev()}", y=0.99)
        else:
            # like phitop, SINGLE plot
            ax = fg.gca()
            f = interp.interp2d(grid["x3"][inds3], grid["x2"][inds2], parm, bounds_error=False)
            parmp = f(yp, xp)
            hi = ax.pcolormesh(xp / 1e3, yp / 1e3, parmp, cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_xlabel('eastward dist. (km)')
            ax.set_ylabel('northward dist. (km)')
            fg.colorbar(hi, ax=ax, label=CB_LBL[name])
            return

        # %% CONVERT TO DISTANCE UP, EAST, NORTH (left panel)
        # JUST PICK AN X3 LOCATION FOR THE MERIDIONAL SLICE PLOT,
        # AND AN ALTITUDE FOR THE LAT./LON. SLICE
        ix3 = lx3 // 2 - 1  # arbitrary slice, to match Matlab
        f = interp.interp2d(grid["x2"][inds2], grid["x1"][inds1], parm[:, :, ix3], bounds_error=False)
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
            xi=np.column_stack((Z3.ravel(), Y3.ravel(), X3.ravel())),
            bounds_error=False,
        ).reshape((1, lyp, lxp))

        parmp = parmp[:, iy, :]  # must be indexed in two steps
        plot23(xp[ix], yp[iy], parmp[0, :, ix], name, cmap, vmin, vmax, fg, axs[1])
        # %% ALT/LAT SLICE (right panel)
        ix2 = lx2 // 2 - 1  # arbitrary slice, to match Matlab
        f = interp.interp2d(grid["x3"][inds3], grid["x1"][inds1], parm[:, ix2, :], bounds_error=False)
        hi = plot13(yp[iy], zp, f(yp, zp)[:, iy], name, cmap, vmin, vmax, fg, axs[2])

        fg.colorbar(hi, ax=axs, aspect=60, pad=0.01)
