#!/usr/bin/env python
"""
plot electric field input to simulation "Efield_inputs" for a single file
"""
import argparse
from matplotlib.pyplot import figure, show
import gemini3d
import numpy as np


def plotVmaxx1it(V: np.ndarray):

    V = V.squeeze()
    fg = figure()
    ax = fg.gca()
    ax.set_title("Vmaxx1it: Potential")
    if V.ndim == 1:
        ax.plot(dat["mlat"], V)
        ax.set_ylabel("Potential [V]")
        ax.set_xlabel("mag. latitude [deg.]")
    elif V.ndim == 2:
        hi = ax.pcolormesh(dat["mlon"], dat["mlat"], V, cmap="bwr")
        ax.set_xlabel("mag. longitude [deg.]")
        ax.set_ylabel("mag. latitude [deg.]")
        fg.colorbar(hi, ax=ax).set_label("potential [V]")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fn", help=".dat or .h5 filename to load directly")
    P = p.parse_args()

    dat = gemini3d.read_Efield(P.fn)

    plotVmaxx1it(dat["Vmaxx1it"][1])

    show()
