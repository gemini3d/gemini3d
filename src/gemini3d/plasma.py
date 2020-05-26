"""
plasma functions
"""
import typing as T
import numpy as np
from pathlib import Path
import shutil
import subprocess
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d, interp2d, interpn

from .config import read_config
from .readdata import readgrid, loadframe
from .base import write_grid, write_state

DictArray = T.Dict[str, T.Any]


def equilibrium_resample(p: T.Dict[str, T.Any], xg: T.Dict[str, T.Any]):
    """
    read and interpolate equilibrium simulation data, writing new
    interpolated grid.
    """

    # %% READ Equilibrium SIMULATION INFO
    peq = read_config(p["eqdir"])
    xgin = readgrid(p["eqdir"])
    # %% END FRAME time of equilibrium simulation
    # this will be the starting time of the new simulation
    t_eq_end = peq["t0"] + peq["tdur"]

    # %% LOAD THE last equilibrium frame
    dat = loadframe(p["eqdir"], t_eq_end)

    # %% sanity check equilibrium simulation input to interpolation
    check_density(dat["ns"][1])
    check_drift(dat["vs"][1])
    check_temperature(dat["Ts"][1])

    # %% DO THE INTERPOLATION
    nsi, vs1i, Tsi = model_resample(xgin, dat["ns"][1], dat["vs"][1], dat["Ts"][1], xg)

    # %% sanity check interpolated variables
    check_density(nsi)
    check_drift(vs1i)
    check_temperature(Tsi)

    # %% WRITE OUT THE GRID
    write_grid(p, xg)

    write_state(t_eq_end, nsi, vs1i, Tsi, p["indat_file"])


def model_resample(
    xgin: DictArray, ns: np.ndarray, vs: np.ndarray, Ts: np.ndarray, xg: DictArray
) -> T.Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ resample a grid
    usually used to upsample an equilibrium simulation grid

    Parameters
    ----------

    xgin: dict
        original grid (usually equilibrium sim grid)
    ns: dict
        number density of species(4D)
    vs: dict
        velocity (4D)
    Ts: dict
        temperature of species (4D)

    Returns
    -------


    """
    # %% NEW GRID SIZES
    lx1, lx2, lx3 = xg["lx"]
    lsp = ns.shape[0]

    # %% ALLOCATIONS
    nsi = np.empty((lsp, lx1, lx2, lx3), dtype=np.float32)
    vsi = np.empty_like(nsi)
    Tsi = np.empty_like(nsi)

    # %% INTERPOLATE ONTO NEWER GRID
    # to avoid IEEE754 rounding issues leading to bounds error,
    # cast the arrays to the same precision,
    # preferring float32 to save disk space and IO time
    X2 = xgin["x2"][2:-2].astype(np.float32)
    X1 = xgin["x1"][2:-2].astype(np.float32)
    X3 = xgin["x3"][2:-2].astype(np.float32)
    x1i = xg["x1"][2:-2].astype(np.float32)
    x2i = xg["x2"][2:-2].astype(np.float32)
    x3i = xg["x3"][2:-2].astype(np.float32)

    if lx3 > 1 and lx2 > 1:
        # 3-D
        print("interpolating grid for 3-D simulation")
        # X2, X1, X3 = np.meshgrid(xgin['x2'][2:-2], xgin['x1'][2:-2], xgin['x3'][2:-2])
        X2i, X1i, X3i = np.meshgrid(x2i, x1i, x3i)
        assert X2i.shape == tuple(xg["lx"])

        for i in range(lsp):
            nsi[i, :, :, :] = interpn(
                points=(X1, X2, X3), values=ns[i, :, :, :], xi=(X1i, X2i, X3i), bounds_error=True
            )
            vsi[i, :, :, :] = interpn(
                points=(X1, X2, X3), values=vs[i, :, :, :], xi=(X1i, X2i, X3i), bounds_error=True
            )
            Tsi[i, :, :, :] = interpn(
                points=(X1, X2, X3), values=Ts[i, :, :, :], xi=(X1i, X2i, X3i), bounds_error=True
            )
    elif lx3 == 1:
        # 2-D east-west
        print("interpolating grid for 2-D simulation in x1, x2")
        # [X2,X1]=meshgrid(xgin.x2(3:end-2),xgin.x1(3:end-2));
        # [X2i,X1i]=meshgrid(xg.x2(3:end-2),xg.x1(3:end-2));
        for i in range(lsp):
            f = interp2d(X2, X1, ns[i, :, :, :], bounds_error=True)
            nsi[i, :, :, :] = f(x2i, x1i)[:, :, None]

            f = interp2d(X2, X1, vs[i, :, :, :], bounds_error=True)
            vsi[i, :, :, :] = f(x2i, x1i)[:, :, None]

            f = interp2d(X2, X1, Ts[i, :, :, :], bounds_error=True)
            Tsi[i, :, :, :] = f(x2i, x1i)[:, :, None]
    elif lx2 == 1:
        # 2-D north-south
        print("interpolating grid for 2-D simulation in x1, x3")
        # original grid, a priori the first 2 and last 2 values are ghost cells
        # on each axis
        #
        # Detect old non-padded grid and workaround
        if np.isclose(xgin["x3"][0], xg["x3"][2], atol=1):
            # old sim, no external ghost cells.
            # Instead of discarding good cells,keep them and say there are
            # new ghost cells outside the grid
            X3 = np.linspace(xgin["x3"][0], xgin["x3"][-1], xgin["lx"][2])
        else:
            # new sim, external ghost cells
            X3 = xgin["x3"][2:-2]

        X1 = xgin["x1"][2:-2]
        # new grid
        x3i = xg["x3"][2:-2].astype(np.float32)
        x1i = xg["x1"][2:-2].astype(np.float32)

        # for each species
        for i in range(lsp):
            f = interp2d(X3, X1, ns[i, :, :, :], bounds_error=True)
            nsi[i, :, :, :] = f(x3i, x1i)[:, None, :]

            f = interp2d(X3, X1, vs[i, :, :, :], bounds_error=True)
            vsi[i, :, :, :] = f(x3i, x1i)[:, None, :]

            f = interp2d(X3, X1, Ts[i, :, :, :], bounds_error=True)
            Tsi[i, :, :, :] = f(x3i, x1i)[:, None, :]

    else:
        raise ValueError("Not sure if this is 2-D or 3-D simulation")

    return nsi, vsi, Tsi


def check_density(n: np.ndarray):

    if not np.isfinite(n).all():
        raise ValueError("non-finite density")
    if (n < 0).any():
        raise ValueError("negative density")
    if n.max() < 1e6:
        raise ValueError("too small maximum density")


def check_drift(v: np.ndarray):

    if not np.isfinite(v).all():
        raise ValueError("non-finite drift")
    if (abs(v) > 10e3).any():
        raise ValueError("excessive drift velocity")


def check_temperature(T: np.ndarray):

    if not np.isfinite(T).all():
        raise ValueError("non-finite temperature")
    if (T < 0).any():
        raise ValueError("negative temperature")
    if T.max() < 500:
        raise ValueError("too cold maximum temperature")


def equilibrium_state(
    p: T.Dict[str, T.Any], xg: DictArray
) -> T.Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    generate (arbitrary) initial conditions for a grid.
    NOTE: only works on symmmetric closed grids!

    [f107a, f107, ap] = activ
    """

    # %% MAKE UP SOME INITIAL CONDITIONS FOR FORTRAN CODE
    mindens = 1e-100

    # %% SLICE THE FIELD IN HALF IF WE ARE CLOSED
    natm = msis_setup(p, xg)

    closeddip = abs(xg["r"][0, 0, 0] - xg["r"][-1, 0, 0]) < 50e3
    # logical flag marking the grid as closed dipole
    if closeddip:
        # closed dipole grid
        #    [~,ialtmax]=max(xg.alt(:,1,1))
        #    lalt=ialtmax
        lalt = xg["lx"][0] // 2
        # FIXME:  needs to work with asymmetric grid...
        alt = xg["alt"][:lalt, :, :]
        lx1 = lalt
        lx2 = xg["lx"][1]
        lx3 = xg["lx"][2]
        Tn = natm[3, :lalt, :, :]
        g = abs(xg["gx1"][:lalt, :, :])
        g = max(g, 1)
        for ix3 in range(lx3):
            for ix2 in range(lx2):
                ialt = abs(g[:, ix2, ix3] - 1).argmin()
                if ialt != lalt:
                    g[ialt:lalt, ix2, ix3] = 1

    else:
        alt = xg["alt"]
        lx1, lx2, lx3 = xg["lx"]
        Tn = natm[3, :, :, :]
        g = abs(xg["gx1"])

    # CONSTANTS
    kb = 1.38e-23
    amu = 1.67e-27

    ns = np.zeros((7, lx1, lx2, lx3))
    for ix3 in range(lx3):
        for ix2 in range(lx2):
            Hf = kb * Tn[:, ix2, ix3] / amu / 16 / g[:, ix2, ix3]
            z0f = 325e3
            He = 2 * kb * Tn[:, ix2, ix3] / amu / 30 / g[:, ix2, ix3]
            z0e = 120e3
            ne = chapmana(alt[:, ix2, ix3], p["nmf"], z0f, Hf) + chapmana(
                alt[:, ix2, ix3], p["nme"], z0e, He
            )
            rho = 1 / 2 * np.tanh((alt[:, ix2, ix3] - 200e3) / 45e3) - 1 / 2 * np.tanh(
                (alt[:, ix2, ix3] - 1000e3) / 200e3
            )

            # has to be .nonzero() as integers not slice is needed.
            inds = (alt[:, ix2, ix3] > z0f).nonzero()[0]
            if len(inds) > 0:
                n0 = p["nmf"]
                #     [n0,ix1]=max(ne);  %in case it isn't exactly z0f
                #     if xg.r(1,1)>xg.r(2,1)
                #         inds=1:ix1;
                #     else
                #         inds=ix1:lx1;
                #     end
                ms = rho[inds] * 16 * amu + (1 - rho[inds]) * amu
                # topside composition only
                H = kb * 2 * Tn[inds, ix2, ix3] / ms / g[inds, ix2, ix3]
                z = alt[inds, ix2, ix3]
                lz = z.size
                iord = np.argsort(z)
                z = z[iord]
                #     z=[z; 2*z(lz)-z(lz-1)];
                z = np.insert(z, 0, z0f)
                integrand = 1 / H[iord]
                integrand = np.append(integrand, integrand[-1])
                #     redheight=intrap(integrand,z);
                redheight = cumtrapz(integrand, z)
                netop = n0 * np.exp(-redheight)
                nesort = np.zeros(lz)
                for iz in range(lz):
                    nesort[iord[iz]] = netop[iz]

                ne[inds] = nesort

            # %% O+
            ns[0, :, ix2, ix3] = rho * ne
            zref = 900e3
            inds0 = alt[:, ix2, ix3] > zref
            if any(inds0):
                iord = np.argsort(alt[:, ix2, ix3])
                altsort = alt[iord, ix2, ix3]
                nsort = ns[0, :, ix2, ix3]
                nsort = nsort[iord]
                #        n0=interpolate(nsort,altsort,zref,'lin','lin');
                f = interp1d(altsort, nsort)
                n0 = f(zref)
                #     [tmp,iref]=min(abs(alt(:,ix2,ix3)-900e3));
                #     if xg.r(1,1)>xg.r(2,1)
                #         inds0=1:iref;
                #     else
                #         inds0=iref:lx1;
                #     end
                #    n0=ns(iref,ix2,ix3,1);
                ms = 16 * amu
                H = kb * 2 * Tn[inds, ix2, ix3] / ms / g[inds, ix2, ix3]
                z = alt[inds0, ix2, ix3]
                lz = z.size
                iord = np.argsort(z)
                z = z[iord]
                #     z=[z; 2*z(lz)-z(lz-1)];
                z = np.insert(z, 0, zref)
                integrand = 1 / H[iord]
                integrand = np.append(integrand, integrand[-1])
                #        redheight=intrap(integrand,z);
                redheight = cumtrapz(integrand, z)
                n1top = n0 * np.exp(-redheight)
                n1sort = np.zeros(lz)
                for iz in range(lz):
                    n1sort[iord[iz]] = n1top[iz]

                ns[0, inds0, ix2, ix3] = n1sort

            # N+
            ns[5, :, ix2, ix3] = 1e-4 * ns[0, :, ix2, ix3]

            inds2 = inds
            inds1 = np.setdiff1d(range(lx1), inds2)

            # MOLECULAR DENSITIES
            nmolc = np.zeros(lx1)
            nmolc[inds1] = (1 - rho[inds1]) * ne[inds1]

            if len(inds2) > 0:
                if xg["r"].ndim == 3:
                    cond = xg["r"][0, 0, 0] > xg["r"][1, 0, 0]
                elif xg["r"].ndim == 2:
                    cond = xg["r"][0, 0] > xg["r"][1, 0]
                else:
                    raise ValueError(
                        "xg['r'] expected to be 3D, possibly with degenerate 2nd or 3rd dimension"
                    )
                if cond:
                    iref = inds1[0]
                else:
                    iref = inds1[-1]

                n0 = nmolc[iref]
                ms = 30.5 * amu
                H = kb * Tn[inds2, ix2, ix3] / ms / g[inds2, ix2, ix3]
                z = alt[inds2, ix2, ix3]
                lz = z.size
                iord = np.argsort(z)
                z = z[iord]
                z = np.append(z, 2 * z[-1] - z[-2])
                integrand = 1 / H[iord]
                integrand = np.append(integrand, integrand[-1])
                #        redheight=intrap(integrand,z);
                redheight = cumtrapz(integrand, z)
                nmolctop = n0 * np.exp(-redheight)
                nmolcsort = np.zeros(lz)
                for iz in range(lz):
                    nmolcsort[iord[iz]] = nmolctop[iz]

                nmolc[inds2] = nmolcsort

            ns[1, :, ix2, ix3] = 1 / 3 * nmolc
            ns[2, :, ix2, ix3] = 1 / 3 * nmolc
            ns[3, :, ix2, ix3] = 1 / 3 * nmolc

            # %% PROTONS
            ns[5, inds2, ix2, ix3] = (1 - rho[inds2]) * ne[inds2]
            z = alt[inds1, ix2, ix3]
            if len(inds2) > 0:
                if cond:
                    iref = inds2[-1]
                else:
                    iref = inds2[0]

                n0 = ns[5, iref, ix2, ix3]
            else:
                iref = alt[:, ix2, ix3].argmax()
                n0 = 1e6

            ns[5, inds1, ix2, ix3] = chapmana(z, n0, alt[iref, ix2, ix3], Hf.mean())

    ns[:6, :, :, :][ns[:6, :, :, :] < mindens] = mindens
    ns[6, :, :, :] = ns[:6, :, :, :].sum(axis=0)

    vsx1 = np.zeros((7, lx1, lx2, lx3))
    Ts = np.tile(Tn[None, :, :, :], [7, 1, 1, 1])

    if closeddip:
        # closed dipole grid
        # FIXME:  This code only works for symmetric grids...
        if 2 * lx1 == xg["lx"][0]:
            ns = np.concatenate((ns, ns[:, ::-1, :, :]), 1)
            Ts = np.concatenate((Ts, Ts[:, ::-1, :, :]), 1)
            vsx1 = np.concatenate((vsx1, vsx1[:, ::-1, :, :]), 1)
        else:
            ns = np.concatenate((ns, ns[:, lx1, :, :], ns[:, ::-1, :, :]), 1)
            Ts = np.concatenate((Ts, Ts[:, lx1, :, :], Ts[:, ::-1, :, :]), 1)
            vsx1 = np.concatenate((vsx1, vsx1[:, lx1, :, :], vsx1[:, ::-1, :, :]), 1)

    return ns, Ts, vsx1


def chapmana(z: np.ndarray, nm: float, z0: float, H: float) -> np.ndarray:
    zref = (z - z0) / H
    ne = nm * np.exp(0.5 * (1 - zref - np.exp(-zref)))

    ne[ne < 1] = 1

    return ne


def msis_setup(p: DictArray, xg: DictArray) -> np.ndarray:
    """calls MSIS Fortran exectuable
    % compiles if not present
    %
    % [f107a, f107, ap] = activ
    %     COLUMNS OF DATA:
    %       1 - ALT
    %       2 - HE NUMBER DENSITY(M-3)
    %       3 - O NUMBER DENSITY(M-3)
    %       4 - N2 NUMBER DENSITY(M-3)
    %       5 - O2 NUMBER DENSITY(M-3)
    %       6 - AR NUMBER DENSITY(M-3)
    %       7 - TOTAL MASS DENSITY(KG/M3)
    %       8 - H NUMBER DENSITY(M-3)
    %       9 - N NUMBER DENSITY(M-3)
    %       10 - Anomalous oxygen NUMBER DENSITY(M-3)
    %       11 - TEMPERATURE AT ALT
    %
    """

    R = Path(__file__).resolve().parent
    builddir = R / "build"
    builddir.mkdir(parents=True, exist_ok=True)
    exe = shutil.which("msis_setup", path=str(builddir))
    if not exe:
        cfg_cmd = ["cmake", "-S", str(R), "-B", str(builddir)]
        print(" ".join(cfg_cmd))
        subprocess.check_call(cfg_cmd)

        build_cmd = ["cmake", "--build", str(builddir), "--target", "msis_setup"]
        print(" ".join(build_cmd))
        subprocess.check_call(build_cmd)

        exe = shutil.which("msis_setup", path=str(builddir))

    if not exe:
        raise FileNotFoundError(f"MSIS setup executable not found in {builddir}")
    # %% SPECIFY SIZES ETC.
    lx1 = xg["lx"][0]
    lx2 = xg["lx"][1]
    lx3 = xg["lx"][2]
    alt = xg["alt"] / 1e3
    glat = xg["glat"]
    glon = xg["glon"]
    lz = lx1 * lx2 * lx3
    # % CONVERT DATES/TIMES/INDICES INTO MSIS-FRIENDLY FORMAT
    t0 = p["t0"]
    doy = int(t0.strftime("%j"))
    UTsec0 = t0.hour * 3600 + t0.minute * 60 + t0.second + t0.microsecond / 1e6

    print("MSIS00 using DOY:", doy)
    yearshort = t0.year % 100
    iyd = yearshort * 1000 + doy
    # %% KLUDGE THE BELOW-ZERO ALTITUDES SO THAT THEY DON'T GIVE INF
    alt[alt <= 0] = 1
    # %% CREATE INPUT FILE FOR FORTRAN PROGRAM
    # don't use NamedTemporaryFile because PermissionError on Windows
    # file_in = tempfile.gettempdir() + "/msis_setup_input.dat"

    # with open(file_in, "w") as f:
    #     np.array(iyd).astype(np.int32).tofile(f)
    #     np.array(UTsec0).astype(np.int32).tofile(f)
    #     np.asarray([p["f107a"], p["f107"], p["Ap"], p["Ap"]]).astype(np.float32).tofile(f)
    #     np.array(lz).astype(np.int32).tofile(f)
    #     np.array(glat).astype(np.float32).tofile(f)
    #     np.array(glon).astype(np.float32).tofile(f)
    #     np.array(alt).astype(np.float32).tofile(f)

    invals = (
        f"{iyd}\n{int(UTsec0)}\n{p['f107a']} {p['f107']} {p['Ap']} {p['Ap']}\n{lz}\n"
        + " ".join(map(str, glat.ravel(order="C")))
        + "\n"
        + " ".join(map(str, glon.ravel(order="C")))
        + "\n"
        + " ".join(map(str, alt.ravel(order="C")))
    )
    # %% CALL MSIS
    # the "-" means to use stdin, stdout
    cmd = [exe, "-", "-", str(lz)]
    print(" ".join(cmd))
    ret = subprocess.check_output(cmd, input=invals, universal_newlines=True)

    Nread = lz * 11

    # old code, from before we used stdout
    # fout_size = Path(file_out).stat().st_size
    # if fout_size != Nread * 4:
    #     raise RuntimeError(f"expected {file_out} size {Nread*4} but got {fout_size}")

    msisdat = np.fromstring(ret, np.float32, Nread, sep=" ").reshape((11, lz), order="F")

    # %% ORGANIZE
    # altitude is a useful sanity check as it's very regular and obvious.
    alt_km = msisdat[0, :].reshape((lx1, lx2, lx3))
    if not np.allclose(alt_km, alt, atol=0.02):  # atol due to precision of stdout ~0.01 km
        raise ValueError("was msis_driver output parsed correctly?")

    nO = msisdat[2, :].reshape((lx1, lx2, lx3))
    nN2 = msisdat[3, :].reshape((lx1, lx2, lx3))
    nO2 = msisdat[4, :].reshape((lx1, lx2, lx3))
    Tn = msisdat[10, :].reshape((lx1, lx2, lx3))
    nN = msisdat[8, :].reshape((lx1, lx2, lx3))

    nNO = 0.4 * np.exp(-3700 / Tn) * nO2 + 5e-7 * nO
    # Mitra, 1968
    nH = msisdat[7, :].reshape((lx1, lx2, lx3))
    natm = np.stack((nO, nN2, nO2, Tn, nN, nNO, nH), 0)

    return natm
