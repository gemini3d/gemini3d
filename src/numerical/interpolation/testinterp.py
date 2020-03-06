#!/usr/bin/env python3
import h5py
from pathlib import Path
import argparse
import sys


def compare_interp(fn: Path, doplot: bool = False):
    fn = Path(fn).expanduser()
    if not fn.is_file():
        print(fn, "not found", file=sys.stderr)
        raise SystemExit(77)

    with h5py.File(fn, "r") as f:
        lx1 = f["lx1"][()]
        lx2 = f["lx2"][()]
        x1 = f["x1"][:]
        x2 = f["x2"][:]

        fx1x2 = f["f"][:]
    print(fx1x2.shape, (lx2, lx1))
    assert fx1x2.shape == (1000, 500) == (lx2, lx1) == (x2.shape, x1.shape)

    if not doplot:
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


# def compare_interp_raw(fn: Path, realbits: int, doplot: bool = False):
#     fn = Path(fn).expanduser()
#     if not fn.is_file():
#         print(fn, "not found", file=sys.stderr)
#         raise SystemExit(77)

#     if realbits == 64:
#         freal = np.float64
#     elif realbits == 32:
#         freal = np.float32
#     else:
#         raise ValueError(f"Unknown realbits {realbits}")

#     with fn.open("r") as f:
#         lx1 = np.fromfile(f, np.int32, 1)[0]
#         lx2 = np.fromfile(f, np.int32, 1)[0]
#         x1 = np.fromfile(f, freal, lx1)
#         x2 = np.fromfile(f, freal, lx2)

#         fx1x2 = np.fromfile(f, freal, lx1 * lx2).reshape((lx1, lx2))
#     assert fx1x2.shape == (500, 1000)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("file")
    p.add_argument("-p", "--plot", help="make plots", action="store_true")
    P = p.parse_args()

    if P.plot:
        from matplotlib.pyplot import figure, show

    compare_interp(P.file, P.plot)

    if P.plot:
        show()
