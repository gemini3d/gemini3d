#!/usr/bin/env python
"""
convert grid .dat to hdf5
"""
import h5py
from numpy import float64, float32, prod, fromfile
from pathlib import Path
import argparse

import gemini3d

read = fromfile
CLVL = 6


def grid2hdf(infn: Path):
    infn = Path(infn).expanduser()
    outfn = infn.with_suffix(".h5")

    lxs = gemini3d.get_simsize(infn.parent / "simsize.dat")
    lgridghost = (lxs[0] + 4) * (lxs[1] + 4) * (lxs[2] + 4)
    gridsizeghost = [lxs[0] + 4, lxs[1] + 4, lxs[2] + 4]

    with infn.open("r") as f, h5py.File(outfn, "w") as h:
        h["lxs"] = lxs
        for i in (1, 2, 3):
            h[f"x{i}"] = read(f, float64, lxs[i - 1] + 4).astype(float32)
            h[f"x{i}i"] = read(f, float64, lxs[i - 1] + 1).astype(float32)
            h[f"dx{i}b"] = read(f, float64, lxs[i - 1] + 3).astype(float32)
            h[f"dx{i}h"] = read(f, float64, lxs[i - 1]).astype(float32)
        for i in (1, 2, 3):
            h.create_dataset(
                f"h{i}",
                data=read(f, float64, lgridghost).reshape(gridsizeghost).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        L = [lxs[0] + 1, lxs[1], lxs[2]]
        for i in (1, 2, 3):
            h.create_dataset(
                f"h{i}x1i",
                data=read(f, float64, prod(L)).reshape(L).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        L = [lxs[0], lxs[1] + 1, lxs[2]]
        for i in (1, 2, 3):
            h.create_dataset(
                f"h{i}x2i",
                data=read(f, float64, prod(L)).reshape(L).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        L = [lxs[0], lxs[1], lxs[2] + 1]
        for i in (1, 2, 3):
            h.create_dataset(
                f"h{i}x3i",
                data=read(f, float64, prod(L)).reshape(L).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        for i in (1, 2, 3):
            h.create_dataset(
                f"gx{i}",
                data=read(f, float64, prod(lxs)).reshape(lxs).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        for k in ("alt", "glat", "glon", "Bmag"):
            h.create_dataset(
                k,
                data=read(f, float64, prod(lxs)).reshape(lxs).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        h.create_dataset(
            "Bincl",
            data=read(f, float64, lxs[1] * lxs[2]).reshape(lxs[1:]).astype(float32),
            compression="gzip",
            compression_opts=CLVL,
        )
        h.create_dataset(
            "nullpts",
            data=read(f, float64, prod(lxs)).reshape(lxs).astype(float32),
            compression="gzip",
            compression_opts=CLVL,
        )
        if f.tell() == infn.stat().st_size:  # not EOF
            return None

        L = [lxs[0], lxs[1], lxs[2], 3]
        for i in (1, 2, 3):
            h.create_dataset(
                f"e{i}",
                data=read(f, float64, prod(L)).reshape(L).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        for k in ("er", "etheta", "ephi"):
            h.create_dataset(
                k,
                data=read(f, float64, prod(L)).reshape(L).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        for k in ("r", "theta", "phi"):
            h.create_dataset(
                k,
                data=read(f, float64, prod(lxs)).reshape(lxs).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )
        if f.tell() == infn.stat().st_size:  # not EOF
            return None

        for k in ("x", "y", "z"):
            h.create_dataset(
                k,
                data=read(f, float64, prod(lxs)).reshape(lxs).astype(float32),
                compression="gzip",
                compression_opts=CLVL,
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("infn", help="Gemini .dat file to convert")
    P = p.parse_args()

    grid2hdf(P.infn)
