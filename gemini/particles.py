import typing as T
import numpy as np
from datetime import timedelta

from .readdata import readgrid


def particles_BCs(p: T.Dict[str, T.Any]):
    """ write particle precipitation to disk """

    # %% GRID ALREADY NEEDS TO BE MADE
    xg = readgrid(p["out_dir"])
    # %% CREATE PRECIPITATION CHARACTERISTICS data
    # number of grid cells.
    # This will be interpolated to grid, so 100x100 is arbitrary
    llon = 100
    llat = 100

    if xg["lx"][1] == 1:  # cartesian
        llon = 1
    elif xg["lx"][2] == 1:
        llat = 1

    # %% TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
    # dtprec is set in config.nml
    time = [p["t0"] + timedelta(seconds=i * p["dtE0"]) for i in range(0, p["tdur"] + p["dtprec"], p["dtprec"])]
    Nt = len(time)

    # %% CREATE PRECIPITATION INPUT DATA
    # Qit: energy flux [mW m^-2]
    # E0it: characteristic energy [eV]

    # did user specify on/off time? if not, assume always on.
    i_on = min(abs(time - p["precip_startsec"])) if "precip_startsec" in p else 0

    i_off = min(abs(time - p["precip_endsec"])) if "precip_endsec" in p else slice(-1)

    pg = precip_grid(xg, p, llat, llon)
    pg["time"] = time

    pg["Q"] = np.empty((Nt, llon, llat))
    pg["E0"] = np.empty((Nt, llon, llat))

    for i in range(i_on, i_off):
        pg["Q"][i, :, :] = precip_gaussian2d(pg)
        pg["E0"][i, :, :] = 5e3

    # %% CONVERT THE ENERGY TO EV
    # E0it = max(E0it,0.100);
    # E0it = E0it*1e3;

    # %% SAVE to files
    # LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
    # FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.
    # THE EFIELD DATA DO NOT NEED TO BE SMOOTHED.

    out_dir = p["out_dir"] / "inputs/prec_inputs/"
    out_dir.mkdir(parents=True, exist_ok=True)


def precip_grid(xg, p, llat: int, llon: int) -> T.Dict[str, T.Any]:

    thetamin = xg["theta"].min()
    thetamax = xg["theta"].max()
    mlatmin = 90 - np.degrees(thetamax)
    mlatmax = 90 - np.degrees(thetamin)
    mlonmin = np.degrees(xg["phi"].min())
    mlonmax = np.degrees(xg["phi"].max())

    # add a 1% buff
    latbuf = 0.01 * (mlatmax - mlatmin)
    lonbuf = 0.01 * (mlonmax - mlonmin)

    E: T.Dict[str, T.Any] = {}
    E["mlat"] = np.linspace(mlatmin - latbuf, mlatmax + latbuf, E["llat"])
    E["mlon"] = np.linspace(mlonmin - lonbuf, mlonmax + lonbuf, E["llon"])
    E["MLON"], E["MLAT"] = np.meshgrid(E["mlon"], E["mlat"])
    # mlonmean = E["mlon"].mean()
    # mlatmean = E["mlat"].mean()

    # %% disturbance width
    mlat_sigma = p["precip_latwidth"] * (mlatmax - mlatmin)
    # to avoid divide by zero below
    E["mlat_sigma"] = max(mlat_sigma, 0.01)
    E["mlon_sigma"] = p["precip_lonwidth"] * (mlonmax - mlonmin)

    return E


def precip_gaussian2d(pg: T.Dict[str, T.Any]) -> np.ndarray:
    return (
        10
        * np.exp(-((pg["MLON"] - pg["mlon_mean"]) ** 2) / (2 * pg["mlon_sigma"] ** 2))
        * np.exp(-((pg["MLAT"] - pg["mlat_mean"]) ** 2) / (2 * pg["mlat_sigma"] ** 2))
    )
