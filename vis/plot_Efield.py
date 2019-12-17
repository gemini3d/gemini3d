#!/usr/bin/env python
"""
plot electric field input to simulation "Efield_inputs"

"""
import argparse
from matplotlib.pyplot import figure, show
import gemini


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fn", help=".dat filename to load directly")
    P = p.parse_args()

    dat = gemini.read_Efield(P.fn)

    fg = figure()
    ax = fg.gca()
    hi = ax.plot(dat["mlat"], dat["Vmaxx1it"].squeeze())
    ax.set_ylabel("Potential [V]")
    ax.set_xlabel("mag. latitude [deg.]")
    show()
