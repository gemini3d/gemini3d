import typing as T
import numpy as np

from .base import write_precip


def particles_BCs(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    """ write particle precipitation to disk """

    # %% CREATE PRECIPITATION CHARACTERISTICS data
    # number of grid cells.
    # This will be interpolated to grid, so 100x100 is arbitrary

    pg: T.Dict[str, T.Any] = {}

    pg["llon"] = 100
    pg["llat"] = 100
    # NOTE: cartesian-specific code
    for k in ("lx", "lxs"):
        if k in xg:
            _, lx2, lx3 = xg[k]
            break
    if lx2 == 1:
        pg["llon"] = 1
    elif lx3 == 1:
        pg["llat"] = 1

    pg = precip_grid(xg, p, pg)

    # %% TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
    # dtprec is set in config.nml
    Nt = (p["tdur"] + p["dtE0"]) // p["dtE0"]
    pg["time"] = [p["t0"] + i * p["dtE0"] for i in range(Nt)]

    # %% CREATE PRECIPITATION INPUT DATA
    # Qit: energy flux [mW m^-2]
    # E0it: characteristic energy [eV]

    # did user specify on/off time? if not, assume always on.
    i_on = min(abs(pg["time"] - p["precip_startsec"])) if "precip_startsec" in p else 0

    i_off = min(abs(pg["time"] - p["precip_endsec"])) if "precip_endsec" in p else Nt

    pg["Q"] = np.empty((Nt, pg["llon"], pg["llat"]))
    pg["E0"] = np.empty((Nt, pg["llon"], pg["llat"]))

    # NOTE: in future, E0 could be made time-dependent in config.nml as 1D array
    for i in range(i_on, i_off):
        pg["Q"][i, :, :] = precip_gaussian2d(pg, p["Qprecip"], p["Qprecip_background"])
        pg["E0"][i, :, :] = p["E0precip"]

    # %% CONVERT THE ENERGY TO EV
    # E0it = max(E0it,0.100);
    # E0it = E0it*1e3;

    # %% SAVE to files
    # LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
    # FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.
    # THE EFIELD DATA DO NOT NEED TO BE SMOOTHED.

    write_precip(pg, p["precdir"], p["format"])


def precip_grid(xg: dict, p: dict, pg: dict) -> T.Dict[str, T.Any]:

    thetamin = xg["theta"].min()
    thetamax = xg["theta"].max()
    mlatmin = 90 - np.degrees(thetamax)
    mlatmax = 90 - np.degrees(thetamin)
    mlonmin = np.degrees(xg["phi"].min())
    mlonmax = np.degrees(xg["phi"].max())

    # add a 1% buff
    latbuf = 0.01 * (mlatmax - mlatmin)
    lonbuf = 0.01 * (mlonmax - mlonmin)

    pg["mlat"] = np.linspace(mlatmin - latbuf, mlatmax + latbuf, pg["llat"])
    pg["mlon"] = np.linspace(mlonmin - lonbuf, mlonmax + lonbuf, pg["llon"])
    # pg["MLON"], pg["MLAT"] = np.meshgrid(pg["mlon"], pg["mlat"])

    # %% disturbance extents
    # avoid divide by zero
    if "precip_latwidth" in p:
        pg["mlat_sigma"] = max(p["precip_latwidth"] * (mlatmax - mlatmin), 0.01)
    if "precip_lonwidth" in p:
        pg["mlon_sigma"] = p["precip_lonwidth"] * (mlonmax - mlonmin)

    return pg


def precip_gaussian2d(pg: T.Dict[str, T.Any], Qpeak: float, Qbackground: float) -> np.ndarray:

    if "mlon_sigma" in pg and "mlat_sigma" in pg:
        Q = (
            Qpeak
            * np.exp(
                -((pg["mlon"][:, None] - pg["mlon"].mean()) ** 2) / (2 * pg["mlon_sigma"] ** 2)
            )
            * np.exp(
                -((pg["mlat"][None, :] - pg["mlat"].mean()) ** 2) / (2 * pg["mlat_sigma"] ** 2)
            )
        )
    elif "mlon_sigma" in pg:
        Q = Qpeak * np.exp(
            -((pg["mlon"][:, None] - pg["mlon"].mean()) ** 2) / (2 * pg["mlon_sigma"] ** 2)
        )
    elif "mlat_sigma" in pg:
        Q = Qpeak * np.exp(
            -((pg["mlat"][None, :] - pg["mlat"].mean()) ** 2) / (2 * pg["mlat_sigma"] ** 2)
        )
    else:
        raise LookupError("precipation must be defined in latitude, longitude or both")

    Q[Q < Qbackground] = Qbackground

    return Q
