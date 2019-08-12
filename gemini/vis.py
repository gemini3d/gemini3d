from datetime import datetime
import numpy as np
import math
from pathlib import Path
import scipy.interpolate as interp
import matplotlib as mpl
from matplotlib.pyplot import figure
import typing

if typing.TYPE_CHECKING:
    import matplotlib.axes as mpla
    import matplotlib.figure as mplf

mpl.rcParams['axes.formatter.limits'] = (-3, 4)
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['axes.formatter.min_exponent'] = 4

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
}


def plotframe(
    time: datetime,
    grid: typing.Dict[str, np.ndarray],
    dat: typing.Dict[str, typing.Any],
    save_dir: Path = None,
    save_ext: str = ".png",
    fg: "mplf.Figure" = None,
):
    """
    if save_dir, plots will not be visible while generating to speed plot writing
    """
    plotfun = grid2plotfun(grid)

    for k in ("ne", "v1", "Ti", "Te", "J1", "v2", "v3", "J2", "J3"):
        if save_dir is None:
            fg = figure(num=k, constrained_layout=True)
        fg.clf()
        plotfun(time, grid, dat[k][1].squeeze(), k, fg)
        if save_dir:
            fg.savefig(save_dir / f"{k}-{time.isoformat().replace(':','')}.{save_ext}")


def grid2plotfun(grid: typing.Dict[str, np.ndarray]):

    minh1 = grid["h1"].min()
    maxh1 = grid["h1"].min()
    if (abs(minh1 - 1) > 1e-4) or (abs(maxh1 - 1) > 1e-4):  # curvilinear grid
        if (grid["lx"][1] > 1) and (grid["lx"][2] > 1):
            plotfun = plot3D_curv_frames_long
        else:
            plotfun = plot2D_curv

    else:  # cartesian grid
        if (grid["lx"][1] > 1) and (grid["lx"][2] > 1):
            plotfun = plot3D_cart_frames_long_ENU
        else:
            plotfun = plot2D_cart

    return plotfun


def plot3D_curv_frames_long(
    time: datetime,
    grid: typing.Dict[str, np.ndarray],
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
):
    pass


def plot2D_curv(
    time: datetime,
    grid: typing.Dict[str, np.ndarray],
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
):
    pass


def plot3D_cart_frames_long_ENU(
    time: datetime,
    grid: typing.Dict[str, np.ndarray],
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
):

    plot_interp(time, grid, parm, name, fg)


def plot2D_cart(
    time: datetime,
    grid: typing.Dict[str, np.ndarray],
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
):

    plot_interp(time, grid, parm, name, fg)


def plot12(
    x: np.ndarray,
    z: np.ndarray,
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
    ax: "mpla.Axes" = None,
):
    if parm.ndim != 2:
        raise ValueError(f"data must have 2 dimensions, you have {parm.shape}")
    if ax is None:
        ax = fg.gca()
    hi = ax.pcolormesh(x / 1e3, z, parm)

    ax.axhline(REF_ALT, color="w", linestyle="--", linewidth=2)
    fg.colorbar(hi, ax=ax, label=CB_LBL[name])

    ax.set_xlabel("eastward dist. (km)")
    ax.set_ylabel("altitude (km)")


def plot13(
    y: np.ndarray,
    z: np.ndarray,
    parm,
    name: str,
    fg: "mplf.Figure",
    ax: "mpla.Axes" = None,
):
    if parm.ndim != 2:
        raise ValueError(f"data must have 2 dimensions, you have {parm.shape}")
    if ax is None:
        ax = fg.gca()
    hi = ax.pcolormesh(y / 1e3, z, parm)

    fg.colorbar(hi, ax=ax, label=CB_LBL[name])

    ax.set_xlabel("northward dist. (km)")
    ax.set_ylabel("altitude (km)")


def plot_interp(
    time: datetime,
    grid: typing.Dict[str, np.ndarray],
    parm: np.ndarray,
    name: str,
    fg: "mplf.Figure",
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
    # %% SIZE OF SIMULATION
    lx1, lx2, lx3 = grid["lx"]
    inds1 = slice(2, lx1 + 2)
    inds2 = slice(2, lx2 + 2)
    inds3 = slice(2, lx3 + 2)
    # %% SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
    meantheta = grid["theta"].mean()
    y = -1 * (grid["theta"] - meantheta)
    x = grid["x2"][inds2] / R_EARTH / math.sin(meantheta)
    z = grid["alt"] / 1e3
    lxp = 500
    lyp = 500
    lzp = 500

    xp = np.linspace(x.min(), x.max(), lxp)
    yp = np.linspace(y.min(), y.max(), lyp)
    zp = np.linspace(z.min(), z.max(), lzp)
    x1plot = zp * 1e3  # upward distance
    x2plot = xp * R_EARTH * math.sin(meantheta)  # eastward distance
    x3plot = yp * R_EARTH  # northward distance
    # %% INTERPOLATE ONTO PLOTTING GRID
    if grid["lx"][2] == 1:  # alt./lon. slice
        ax = fg.gca()
        ax.set_title(time.isoformat())
        # meridional meshgrid, this defines the grid for plotting
        # slice expects the first dim. to be "y" ("z" in the 2D case)
        fmp = interp.interp2d(grid["x2"][inds2], grid["x1"][inds1], parm)
        parmp = fmp(x2plot, x1plot)
        # %% CONVERT ANGULAR COORDINATES TO MLAT,MLON
        xp *= R_EARTH * math.sin(meantheta)  # eastward ground distance
        i = np.argsort(xp)
        plot12(xp[i], zp, parmp[:, i], name, fg)
    elif grid["lx"][1] == 1:  # alt./lat. slice
        ax = fg.gca()
        ax.set_title(time.isoformat())
        # so north dist, east dist., alt.
        # slice expects the first dim. to be "y"
        fmp = interp.interp2d(grid["x3"][inds3], grid["x1"][inds1], parm)
        parmp = fmp(x3plot, x1plot).reshape((lzp, lyp))
        # %% CONVERT ANGULAR COORDINATES TO MLAT,MLON
        i = np.argsort(yp)
        plot13(yp[i] * R_EARTH, zp, parmp[:, i], name, fg)
    else:  # 3-panel plot, vs. single-panel plots of 2-D cases
        fg.set_size_inches((18, 5))
        axs = fg.subplots(1, 3)
        fg.suptitle(time.isoformat(), y=0.99)

        xp *= 1e6
        yp *= 1e6
        # %% CONVERT TO DISTANCE UP, EAST, NORTH
        # JUST PICK AN X3 LOCATION FOR THE MERIDIONAL SLICE PLOT,
        # AND AN ALTITUDE FOR THE LAT./LON. SLICE
        fmp = interp.interp2d(
            grid["x2"][inds2], grid["x1"][inds1], parm[:, :, lx3 // 2]
        )
        # slice expects the first dim. to be "y" ("z" in the 2D case)
        parmp = fmp(x2plot, zp)
        # CONVERT ANGULAR COORDINATES TO MLAT,MLON
        ix = np.argsort(xp)
        iy = np.argsort(yp)
        plot12(xp[ix], zp, parmp[:, ix], name, fg, axs[0])
        # %% LAT./LONG. SLICE COORDINATES
        zp2 = np.array([REF_ALT - 10, REF_ALT, REF_ALT + 10])
        X3, Y3, Z3 = np.meshgrid(x2plot, x3plot, zp2 * 1e3)
        # transpose: so north dist, east dist., alt.
        parmp = interp.interpn(
            (
                grid["x2"][inds2],
                grid["x3"][
                    inds3
                ],  # this is northward distance - again backwards from yp
                grid["x1"][inds1],
            ),
            values=np.transpose(parm, [2, 1, 0]),
            xi=np.column_stack(
                (X3.ravel(), Y3.ravel(), Z3.ravel())
            ),  # slice expects the first dim. to be "y"
            bounds_error=False,
        ).reshape((lyp, lxp, len(zp2)))

        parmp = parmp[iy, :, :]  # must be indexed in two steps
        plot13(xp[ix], yp[iy], parmp[:, ix, 2], name, fg, axs[1])
        # %% ALT/LAT SLICE
        fmp = interp.interp2d(
            grid["x3"][inds3], grid["x1"][inds1], parm[:, lx2 // 2, :]
        )
        parmp = fmp(x3plot, x1plot)
        plot13(yp[iy], zp, parmp[:, iy], name, fg, axs[2])
