#!/usr/bin/env python3
import h5py
from pathlib import Path
import argparse


def compare_interp(fn: Path, doplot: bool = False) -> None:

    with h5py.File(Path(fn).expanduser(), "r") as f:
        lx1 = f["/lx1"][()]
        lx2 = f["/lx2"][()]
        lx3 = f["/lx3"][()]
        x1 = f["/x1"][:]
        x2 = f["/x2"][:]
        x3 = f["/x3"][:]

        fx1x2x3 = f["/f"][:]

    # NOTE: C order on disk
    assert fx1x2x3.shape == (
        x3.size,
        x2.size,
        x1.size,
    ), f"{x3.size} {x2.size} {x1.size} != {fx1x2x3.shape}   {lx3} {lx2} {lx1}"

    if not doplot:
        return None

    f1 = figure()
    ax = f1.gca()
    ix3 = lx3 // 2
    hi = ax.pcolormesh(x2, x1, fx1x2x3[ix3, :, :].transpose(), shading="nearest")
    ax.set_xlabel("$x_2$")
    ax.set_ylabel("$x_1$")
    f1.colorbar(hi, ax=ax).set_label("$f(x_2, x_1)$" + f" {lx2}, {lx1}")
    ax.set_title("$f(x_3, x_1) x_3$=" + f"{ix3}")

    f2 = figure()
    ax = f2.gca()
    ix2 = lx2 // 2
    hi = ax.pcolormesh(x3, x1, fx1x2x3[:, ix2, :].transpose(), shading="nearest")
    ax.set_xlabel("$x_3$")
    ax.set_ylabel("$x_1$")
    f2.colorbar(hi, ax=ax).set_label("$f(x_3, x_1)$" +  f" {lx3}, {lx1}")
    ax.set_title("$f(x_3, x_1) x_2$= " + f"{ix2}")

    f3 = figure()
    ax = f3.gca()
    ix1 = lx1 // 2
    hi = ax.pcolormesh(x2, x3, fx1x2x3[:, :, ix1], shading="nearest")
    ax.set_xlabel("$x_2$")
    ax.set_ylabel("$x_3$")
    f3.colorbar(hi, ax=ax).set_label("$f(x_2, x_3)$" + f" {lx2}, {lx3}")
    ax.set_title("$f(x_2, x_3) x_1$= " + f"{ix1}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("file", help="data file to load")
    p.add_argument("-p", "--plot", help="make plots", action="store_true")
    P = p.parse_args()

    if P.plot:
        from matplotlib.pyplot import figure, show

    compare_interp(P.file, P.plot)

    if P.plot:
        show()

    print("OK: test interp 3d")
