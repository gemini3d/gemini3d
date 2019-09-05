# -*- coding: utf-8 -*-
"""
convert Gemini .dat to .h5
"""
import h5py
from pathlib import Path
import argparse
import typing
from numpy import float32

import gemini

LSP = 7
CLVL = 6


def dat2hdf(infn: Path, outfn: Path):

    P = gemini.readconfig(infn.parent / "inputs/config.ini")
    if P["flagoutput"] == 1:
        convert3d_curv(infn, outfn, P["lxs"])
    elif P["flagoutput"] == 2:
        convert3d_curvavg(infn, outfn, P["lxs"])
    else:
        raise ValueError(f"not sure how to convert {infn}")


def convert3d_curv(infn: Path, outfn: Path, lxs: typing.Sequence[int]):

    with infn.open("r") as f, h5py.File(outfn, "w") as h:
        h["time"] = gemini.read_time(f).isoformat()
        h.create_dataset(
            "ns",
            data=gemini.read4D(f, LSP, lxs).astype(float32),
            chunks=(1, *lxs[1:], LSP),
            compression="gzip",
            compression_opts=CLVL,
        )
        h.create_dataset(
            "vs1",
            data=gemini.read4D(f, LSP, lxs).astype(float32),
            chunks=(1, *lxs[1:], LSP),
            compression="gzip",
            compression_opts=CLVL,
        )
        h.create_dataset(
            "Ts",
            data=gemini.read4D(f, LSP, lxs).astype(float32),
            chunks=(1, *lxs[1:], LSP),
            compression="gzip",
            compression_opts=CLVL,
        )

        for p in ("J1", "J2", "J3", "v2", "v3"):
            h.create_dataset(
                p, data=gemini.read3D(f, lxs).astype(float32), chunks=(1, *lxs[1:]), compression="gzip", compression_opts=CLVL
            )

        h.create_dataset("Phitop", data=gemini.read2D(f, lxs).astype(float32), compression="gzip", compression_opts=CLVL)


def convert3d_curvavg(infn: Path, outfn: Path, lxs: typing.Sequence[int]):

    with infn.open("r") as f, h5py.File(outfn, "w") as h:
        h["time"] = gemini.read_time(f).isoformat()

        for p in ("ne", "v1", "Ti", "Te", "J1", "J2", "J3", "v2", "v3"):
            h.create_dataset(
                p, data=gemini.read3D(f, lxs).astype(float32), chunks=(1, *lxs[1:]), compression="gzip", compression_opts=CLVL
            )

        h.create_dataset("Phitop", data=gemini.read2D(f, lxs).astype(float32), compression="gzip", compression_opts=CLVL)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("indir", help="Gemini .dat file directory")
    p.add_argument("-o", "--outdir", help="directory to write HDF5 files")
    P = p.parse_args()

    indir = Path(P.indir).expanduser()
    outdir = Path(P.outdir).expanduser() if P.outdir else indir

    for infile in indir.glob("*.dat"):
        outfile = outdir / (infile.stem + ".h5")
        print(infile, "=>", outfile)
        dat2hdf(infile, outfile)
