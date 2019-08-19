#!/usr/bin/env python3
import numpy as np
from pathlib import Path
import argparse
import sys

try:
    from matplotlib.pyplot import figure, show
except (ImportError, RuntimeError):
    figure = show = None


def compare_interp(fn: Path):
    fn = Path(fn).expanduser()
    if not fn.is_file():
        print(fn, "not found", file=sys.stderr)
        raise SystemExit(77)

    with fn.open("r") as f:
        lx1 = np.fromfile(f, np.int32, 1)[0]
        lx2 = np.fromfile(f, np.int32, 1)[0]
        x1 = np.fromfile(f, np.float64, lx1)
        x2 = np.fromfile(f, np.float64, lx2)

        fx1x2 = np.fromfile(f, np.float64, lx1 * lx2).reshape((lx1, lx2))
        assert fx1x2.shape == (500, 1000)

    if show is None:
        return None

    fg = figure()
    ax = fg.gca()

    if lx2 == 1:
        ax.plot(x1, fx1x2)
        ax.set_xlabel("x_1")
        ax.set_ylabel("fx1x2")
        ax.set_title("1-D interp")
    else:
        hi = ax.pcolormesh(x2, x1, fx1x2)
        ax.set_xlabel("x_2")
        ax.set_ylabel("x_1")
        fg.colorbar(hi, ax=ax).set_label("fx1x2")
        ax.set_title("2-D interp")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("file")
    P = p.parse_args()

    compare_interp(P.file)
