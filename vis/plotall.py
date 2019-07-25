#!/usr/bin/env python
from argparse import ArgumentParser
import concurrent.futures
from pathlib import Path
from copy import deepcopy
import itertools
import oct2py


def frame(direc: Path, ymd, UTsec, saveplots):
    with oct2py.Oct2Py() as octave:
        octave.plotframe(str(direc), ymd, UTsec, saveplots)


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

    direc = Path(p.direc).expanduser()
    config = direc / "inputs/config.ini"

    # works single-threaded, was just a first test
    #    with oct2py.Oct2Py() as octave:
    #       octave.plotall(p.direc, p.saveplots)
    # %% setup
    with oct2py.Oct2Py() as octave:
        octave.addpath("../script_utils")
        ymd0, UTsec0, tdur, dtout = octave.readconfig(str(config), nout=4)

        Nt = int((UTsec0 + tdur - UTsec0) // dtout) + 1

        ymd = deepcopy(ymd0)
        UTsec = deepcopy(UTsec0)
        # %% setup threading
        ymds = []
        UTsecs = []
        for _ in range(Nt):
            ymd, UTsec = octave.dateinc(dtout, ymd, UTsec, nout=2)
            ymds.append(ymd)
            UTsecs.append(UTsec)
    # %% threading
    # set max_workers as high as your PC can handle
    with concurrent.futures.ThreadPoolExecutor(max_workers=p.max_threads) as exc:
        exc.map(
            frame, itertools.repeat(direc), ymds, UTsecs, itertools.repeat(p.saveplots)
        )


if __name__ == "__main__":
    main()
