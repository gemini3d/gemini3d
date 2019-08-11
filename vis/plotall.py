#!/usr/bin/env python
"""
plots simulation output--a simple example
"""
from argparse import ArgumentParser
from pathlib import Path
import matplotlib.pyplot as mpl

import gemini.readdata as grd
import gemini.vis as vis


def main():
    p = ArgumentParser()
    p.add_argument("direc", help="directory to plot")
    p.add_argument(
        "-s",
        "--saveplots",
        help="save plots: directory, type.  e.g. -s /tmp png  or  -s /tmp eps",
        nargs=2,
    )
    p = p.parse_args()

    direc = Path(p.direc).expanduser().resolve(strict=True)
    if p.saveplots:
        save_dir, save_ext = p.saveplots
        save_dir = Path(save_dir).expanduser()
        save_dir.mkdir(parents=True, exist_ok=True)
        from matplotlib.figure import Figure

        fg = Figure()
    else:
        save_dir = save_ext = None
        fg = None

    params = grd.readconfig(direc / "inputs/config.ini")
    t0 = params["t0"]
    times = grd.datetime_range(t0, t0 + params["tdur"], params["dtout"])

    grid = grd.readgrid(direc / "inputs/simgrid.dat")

    for t in times:
        dat = grd.loadframe(direc, t)
        vis.plotframe(t, grid, dat, save_dir, save_ext, fg)
        if p.saveplots:
            print(f"saving {t} to {save_dir}")
        else:
            mpl.draw()
            mpl.pause(1)


if __name__ == "__main__":
    main()
