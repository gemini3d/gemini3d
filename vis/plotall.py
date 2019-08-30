#!/usr/bin/env python3
"""
plots simulation output--a simple example
"""
from argparse import ArgumentParser
from pathlib import Path
import logging
import matplotlib.pyplot as mpl

import gemini.readdata as grd
import gemini.vis as vis


def main():
    p = ArgumentParser()
    p.add_argument("direc", help="directory to plot")
    p.add_argument("-s", "--saveplots", help="save plots to data directory", action="store_true")
    p = p.parse_args()

    direc = Path(p.direc).expanduser().resolve(strict=True)
    if p.saveplots:
        from matplotlib.figure import Figure

        fg = Figure(tight_layout=True)
        save_dir = direc / "plots"
        save_dir.mkdir(parents=True, exist_ok=True)
    else:
        fg = None
        save_dir = None

    params = grd.readconfig(direc / "inputs/config.ini")
    t0 = params["t0"]
    times = grd.datetime_range(t0, t0 + params["tdur"], params["dtout"])

    grid = grd.readgrid(direc / "inputs/simgrid.dat")

    for t in times:
        try:
            dat = grd.loadframe(direc, t)
        except Exception as err:
            logging.error(f"{t} in {direc} not loadable: {err}")
            continue

        vis.plotframe(t, grid, dat, save_dir, fg)
        if not p.saveplots:
            mpl.draw()
            mpl.pause(1)
        else:
            print(f"saving {t} to {save_dir}")


if __name__ == "__main__":
    main()
