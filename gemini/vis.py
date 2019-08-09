from datetime import datetime
import numpy as np
import math
import scipy.interpolate as interp
from matplotlib.pyplot import figure, draw, pause

R_EARTH = 6370e3
REF_ALT = 300  # km


def plotframe(time: datetime, grid: dict, dat: dict):
    plotfun = grid2plotfun(grid)

    for k in ("ne", "v1", "Ti", "Te", "J1", "v2", "v3", "J2", "J3"):
        fg = figure(k)
        ax = fg.gca()
        plotfun(time, grid, dat[k][1].squeeze(), k, fg, ax)
        draw()
        pause(1)


def grid2plotfun(grid: dict):

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
    time: datetime, grid: dict, parm: dict, name: str, fg=None, ax=None
):
    pass


def plot2D_curv(time: datetime, grid: dict, parm: dict, name: str, fg=None, ax=None):
    pass


def plot3D_cart_frames_long_ENU(
    time: datetime, grid: dict, parm: dict, name: str, fg=None, ax=None
):
    pass


def plot2D_cart(time: datetime, grid: dict, parm: dict, name: str, fg=None, ax=None):

    if fg is None:
        fg = figure()
        ax = fg.gca()
    # %% SIZE OF SIMULATION
    lx1, lx2, lx3 = grid["lx"]
    inds1 = slice(2, lx1 + 2)
    inds2 = slice(2, lx2 + 2)
    inds3 = slice(2, lx3 + 2)
    # %% SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
    meantheta = grid["theta"].mean()
    """
    this is a mag colat. coordinate and is only used for defining
    grid in linspaces below, runs backward from north distance,
    hence the negative sign
    """
    y = -1 * (grid["theta"] - meantheta)
    x = grid["x2"][inds2] / R_EARTH / math.sin(meantheta)
    z = grid["alt"] / 1e3
    lxp = 500
    lyp = 500
    lzp = 500

    xp = np.linspace(x.min(), x.max(), lxp)  # eastward distance (rads.)
    # should be interpreted as northward distance (in rads.).
    # Irrespective of ordering of xg.theta, this will be monotonic increasing!!!
    yp = np.linspace(y.min(), y.max(), lyp)
    zp = np.linspace(z.min(), z.max(), lzp)  # altitude (meters)
    # %% INTERPOLATE ONTO PLOTTING GRID
    if grid["lx"][2] == 1:  # alt./lon. slice
        # meridional meshgrid, this defines the grid for plotting
        x1plot = zp * 1e3  # upward distance
        x2plot = xp * R_EARTH * math.sin(meantheta)  # eastward distance
        # slice expects the first dim. to be "y" ("z" in the 2D case)
        fmp = interp.interp2d(grid["x2"][inds2], grid["x1"][inds1], parm)
        parmp = fmp(x2plot, x1plot).reshape((lzp, lxp))

    elif grid["lx"][1] == 1:  # alt./lat. slice
        x1plot = zp * 1e3  # upward distance
        x3plot = yp * R_EARTH  # northward distance;

        # so north dist, east dist., alt.
        # slice expects the first dim. to be "y"
        fmp3 = interp.interp2d(grid["x3"][inds3], grid["x1"][inds1], parm)
        parmp3 = fmp3(x3plot, x1plot).reshape((lzp, lyp))

    # %% CONVERT ANGULAR COORDINATES TO MLAT,MLON
    if lx2 == 1:
        yp *= R_EARTH
        i = np.argsort(yp)
        yp = yp[i]
        parmp3 = parmp3[:, i]

    if lx3 == 1:
        xp *= R_EARTH * math.sin(meantheta)  # eastward ground distance
        i = np.argsort(xp)
        xp = xp[i]
        parmp = parmp[:, i]

    if lx3 == 1:
        plot12(xp, zp, parmp, fg, ax)
    elif lx2 == 1:
        plot13(yp, zp, parmp3, fg, ax)


def plot12(x, z, parm, fg, ax):
    hi = ax.pcolormesh(x / 1e3, z, parm)

    ax.axhline(REF_ALT, color="w", linestyle="--", linewidth=2)
    fg.colorbar(hi, ax=ax)

    ax.set_xlabel("eastward dist. (km)")
    ax.set_ylabel("altitude (km)")


def plot13(y, z, parm, fg, ax):
    hi = ax.pcolormesh(y / 1e3, z, parm)

    fg.colorbar(hi, ax=ax)

    ax.set_xlabel("northward dist. (km)")
    ax.set_ylabel("altitude (km)")
