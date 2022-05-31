#!/usr/bin/env python3
"""
do a simple test of Fang 2008

This program is typically used from CMake as a unit test for Fang ionization.

It can also be used by a human to plot ionization profiles to compare with the original
Fang 2008 and Fang 2010 papers.
"""
import numpy as np
import subprocess
from pathlib import Path
import shutil
import sys
import io
import argparse

Rb = Path(__file__).resolve().parents[2] / "build/src/ionization"


def checker(exe: str, doplot: bool, params: dict = None):
    if not exe:
        if not Rb.is_dir():
            raise FileNotFoundError(
                f"build directory does not exist, did you build Gemini with CMake?  {Rb}"
            )
        exe = shutil.which("test_fang", path=str(Rb))

    if not shutil.which(exe):
        print("test_fang executable not found:", exe, file=sys.stderr)
        raise SystemExit(77)

    if params:
        ret = subprocess.check_output(
            [
                exe,
                params["Q0"],
                params["E0"],
                params["alt_range"],
                params["f107"],
                params["f107a"],
                params["Ap"],
                params["glat"],
                params["glon"],
                params["doy"],
                params["utsec"],
            ],
            universal_newlines=True,
        )
    else:
        ret = subprocess.check_output(exe, universal_newlines=True)

    keV = list(map(float, ret.split("\n")[0].split()[1:]))
    dat = np.loadtxt(io.StringIO(ret), skiprows=1, max_rows=191)
    alt_km = dat[:, 0]
    ionization_rates08 = dat[:, 1:]

    dat = np.loadtxt(io.StringIO(ret), skiprows=191 + 3, max_rows=191)
    ionization_rates10 = dat[:, 1:]

    if not params:
        # spot check values
        assert np.isclose(ionization_rates08[89, 0], 2214.052, atol=0.001), "E0: 100eV"
        assert np.isclose(ionization_rates08[17, 4], 9579.046, atol=0.001), "E0: 1MeV"

        assert np.isclose(ionization_rates10[89, 0], 1192.002, atol=0.001), "Emono: 100eV"
        assert np.isclose(ionization_rates10[17, 4], 778.655, atol=0.001, rtol=0.001), "Emono: 1MeV"

    if not doplot:
        return

    fg = figure()
    axs = fg.subplots(1, 2, sharey=True)
    fg.suptitle(r"Ap=5 f107=50 Midnight MLT 60$^\circ$ lat.")
    # %% Fang 2008 plot
    ax = axs[0]
    for i, e in enumerate(keV):
        ax.semilogx(ionization_rates08[:, i], alt_km, label=str(e))
    ax.set_ylabel("altitude [km]")
    ax.set_xlabel("Total ionization rate [cm$^{-3}$ s$^{-1}$]")
    ax.grid(True)
    ax.set_title(r"Figure 3 of Fang 2008 by $E_0$ [keV]")
    ax.legend(loc="best")
    ax.set_xlim(10, 1e5)
    # %% Fang 2010 plot
    ax = axs[1]
    for i, e in enumerate(keV):
        ax.semilogx(ionization_rates10[:, i], alt_km, label=str(e))
    ax.set_xlabel("Total ionization rate [cm$^{-3}$ s$^{-1}$]")
    ax.grid(True)
    ax.set_title(r"Figure 2 of Fang 2010 by $E_{mono}$ [keV]")
    ax.legend(loc="best")
    ax.set_xlim(10, 1e5)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("exe", help="path to test_fang executable (as called by CMake)", nargs="?")
    p.add_argument("-p", "--plot", help="make plots", action="store_true")
    P = p.parse_args()

    if P.plot:
        from matplotlib.pyplot import figure, show

    checker(P.exe, P.plot)

    if P.plot:
        show()
