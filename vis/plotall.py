#!/usr/bin/env python
"""
plots simulation output--a simple example
"""
from argparse import ArgumentParser
from pathlib import Path

import gemini.readdata as grd
import gemini.vis as vis


def main():
    p = ArgumentParser()
    p.add_argument("direc", help="directory to plot")
    p.add_argument(
        "max_threads",
        help="maximum number of threads to use (large grids use lots of RAM)",
        nargs="?",
        type=int,
        default=5,
    )
    p.add_argument(
        "-s",
        "--saveplots",
        help="plot type to save (png, eps) [default png]",
        nargs="+",
        default=["png"],
    )
    p = p.parse_args()

    direc = Path(p.direc).expanduser().resolve(strict=True)

    params = grd.readconfig(direc / "inputs/config.ini")
    t0 = params["t0"]
    times = grd.datetime_range(t0, t0 + params["tdur"], params["dtout"])

    grid = grd.readgrid(direc / "inputs/simgrid.dat")

    for t in times:
        dat = grd.loadframe(direc, t)
        vis.plotframe(t, grid, dat)


if __name__ == "__main__":
    main()
