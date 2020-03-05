import typing as T
import numpy as np
import logging
from scipy.special import erf

from .base import write_Efield


def Efield_BCs(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]) -> T.Dict[str, T.Any]:
    """ generate E-field """

    E: T.Dict[str, T.Any] = {}

    E["Efield_outdir"] = p["out_dir"] / "inputs/Efield_inputs"
    E["Efield_outdir"].mkdir(parents=True, exist_ok=True)

    # %% READ IN THE SIMULATION INFO
    for k in ("lx", "lxs"):
        if k in xg:
            lx1, lx2, lx3 = xg[k]
            break

    # %% CREATE ELECTRIC FIELD DATASET
    E["llon"] = 100
    E["llat"] = 100
    # NOTE: cartesian-specific code
    if lx2 == 1:
        E["llon"] = 1
    elif lx3 == 1:
        E["llat"] = 1

    thetamin = xg["theta"].min()
    thetamax = xg["theta"].max()
    mlatmin = 90 - np.degrees(thetamax)
    mlatmax = 90 - np.degrees(thetamin)
    mlonmin = np.degrees(xg["phi"].min())
    mlonmax = np.degrees(xg["phi"].max())

    # add a 1% buff
    latbuf = 0.01 * (mlatmax - mlatmin)
    lonbuf = 0.01 * (mlonmax - mlonmin)
    E["mlat"] = np.linspace(mlatmin - latbuf, mlatmax + latbuf, E["llat"])
    E["mlon"] = np.linspace(mlonmin - lonbuf, mlonmax + lonbuf, E["llon"])
    # E["MLON"], E["MLAT"] = np.meshgrid(E["mlon"], E["mlat"])
    mlonmean = E["mlon"].mean()
    mlatmean = E["mlat"].mean()

    # %% WIDTH OF THE DISTURBANCE
    mlatsig = p["Efield_fracwidth"] * (mlatmax - mlatmin)
    mlonsig = p["Efield_fracwidth"] * (mlonmax - mlonmin)
    sigx2 = p["Efield_fracwidth"] * (xg["x2"].max() - xg["x2"].min())
    sigx3 = p["Efield_fracwidth"] * (xg["x3"].max() - xg["x3"].min())
    # %% TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
    Nt = (p["tdur"] + p["dtE0"]) // p["dtE0"]
    E["time"] = [p["t0"] + i * p["dtE0"] for i in range(Nt)]
    # time given in file is the seconds from UTC midnight

    # %% CREATE DATA FOR BACKGROUND ELECTRIC FIELDS
    # assign to zero in case not specifically assigned
    E["Exit"] = np.zeros((Nt, E["llon"], E["llat"]))
    E["Eyit"] = np.zeros((Nt, E["llon"], E["llat"]))

    # %% CREATE DATA FOR BOUNDARY CONDITIONS FOR POTENTIAL SOLUTION
    p["flagdirich"] = 1
    # if 0 data is interpreted as FAC, else we interpret it as potential
    E["Vminx1it"] = np.zeros((Nt, E["llon"], E["llat"]))
    E["Vmaxx1it"] = np.zeros((Nt, E["llon"], E["llat"]))
    # these are just slices
    E["Vminx2ist"] = np.zeros((Nt, E["llat"]))
    E["Vmaxx2ist"] = np.zeros((Nt, E["llat"]))
    E["Vminx3ist"] = np.zeros((Nt, E["llon"]))
    E["Vmaxx3ist"] = np.zeros((Nt, E["llon"]))

    if p["Etarg"] > 1:
        logging.warning(f"Etarg units V/m -- is {p['Etarg']} V/m realistic?")

    # NOTE: h2, h3 have ghost cells, so we use lx1 instead of -1 to index
    # pk is a scalar.
    if lx3 == 1:
        # east-west
        pk = p["Etarg"] * sigx2 * xg["h2"][lx1, lx2 // 2, 0] * np.sqrt(np.pi) / 2
    elif lx2 == 1:
        # north-south
        pk = p["Etarg"] * sigx3 * xg["h3"][lx1, 0, lx3 // 2] * np.sqrt(np.pi) / 2
    else:
        # 3D
        pk = p["Etarg"] * sigx2 * xg["h2"][lx1, lx2 // 2, 0] * np.sqrt(np.pi) / 2

    assert pk.ndim == 0, "pk is a scalar"

    for i in range(Nt):
        if lx2 == 1:
            E["Vmaxx1it"][i, :, :] = pk * erf((E["mlat"] - mlatmean) / mlatsig)[None, :]
        elif lx3 == 1:
            E["Vmaxx1it"][i, :, :] = pk * erf((E["mlon"] - mlonmean) / mlonsig)[:, None]
        else:
            E["Vmaxx1it"][i, :, :] = pk * (
                erf((E["mlon"] - mlonmean) / mlonsig)[:, None] * erf((E["mlat"] - mlatmean) / mlatsig)[None, :]
            )

    # %% check for NaNs
    # this is also done in Fortran, but just to help ensure results.
    check_finite(E["Exit"], "Exit")
    check_finite(E["Eyit"], "Eyit")
    check_finite(E["Vminx1it"], "Vminx1it")
    check_finite(E["Vmaxx1it"], "Vmaxx1it")
    check_finite(E["Vminx2ist"], "Vminx2ist")
    check_finite(E["Vmaxx2ist"], "Vmaxx2ist")
    check_finite(E["Vminx3ist"], "Vminx3ist")
    check_finite(E["Vmaxx3ist"], "Vmaxx3ist")

    # %% SAVE THESE DATA TO APPROPRIATE FILES
    # LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
    # FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.
    # THE EFIELD DATA DO NOT TYPICALLY NEED TO BE SMOOTHED.
    write_Efield(p, E)

    return E


def check_finite(v: np.ndarray, name: str):

    i = ~np.isfinite(v)
    if i.any():
        raise ValueError(f"{np.count_nonzero(i)} NaN in {name} at {i.nonzero()}")
