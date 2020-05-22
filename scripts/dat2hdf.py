# -*- coding: utf-8 -*-
"""
convert Gemini .dat to .h5
"""
import h5py
from pathlib import Path
import argparse
import typing
from numpy import float32

import gemini3d

LSP = 7
CLVL = 6


def dat2hdf(infn: Path, outfn: Path):

    P = gemini3d.read_config(infn.parent / "inputs")
    if P["flagoutput"] == 1:
        convert3d_curv(infn, outfn, P["lxs"])
    elif P["flagoutput"] == 2:
        convert3d_curvavg(infn, outfn, P["lxs"])
    else:
        raise ValueError(f"not sure how to convert {infn}")


def convert3d_curv(infn: Path, outfn: Path, lxs: typing.Sequence[int]):

    with infn.open("r") as f, h5py.File(outfn, "w") as h:
        h["time"] = gemini3d.read_time(f).isoformat()
        h.create_dataset(
            "ns",
            data=gemini3d.read4D(f, LSP, lxs).astype(float32),
            chunks=(1, *lxs[1:], LSP),
            compression="gzip",
            compression_opts=CLVL,
        )
        h.create_dataset(
            "vs1",
            data=gemini3d.read4D(f, LSP, lxs).astype(float32),
            chunks=(1, *lxs[1:], LSP),
            compression="gzip",
            compression_opts=CLVL,
        )
        h.create_dataset(
            "Ts",
            data=gemini3d.read4D(f, LSP, lxs).astype(float32),
            chunks=(1, *lxs[1:], LSP),
            compression="gzip",
            compression_opts=CLVL,
        )

        for p in ("J1", "J2", "J3", "v2", "v3"):
            h.create_dataset(
                p,
                data=gemini3d.read3D(f, lxs).astype(float32),
                chunks=(1, *lxs[1:]),
                compression="gzip",
                compression_opts=CLVL,
            )

        h.create_dataset(
            "Phitop",
            data=gemini3d.read2D(f, lxs).astype(float32),
            compression="gzip",
            compression_opts=CLVL,
        )


def convert3d_curvavg(infn: Path, outfn: Path, lxs: typing.Sequence[int]):

    with infn.open("r") as f, h5py.File(outfn, "w") as h:
        h["time"] = gemini3d.read_time(f).isoformat()

        for p in ("ne", "v1", "Ti", "Te", "J1", "J2", "J3", "v2", "v3"):
            h.create_dataset(
                p,
                data=gemini3d.read3D(f, lxs).astype(float32),
                chunks=(1, *lxs[1:]),
                compression="gzip",
                compression_opts=CLVL,
            )

        h.create_dataset(
            "Phitop",
            data=gemini3d.read2D(f, lxs).astype(float32),
            compression="gzip",
            compression_opts=CLVL,
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("indir", help="Gemini .dat file directory")
    p.add_argument("-o", "--outdir", help="directory to write HDF5 files")
    p.add_argument(
        "--delete", help="delete original file if conversion successful", action="store_true"
    )
    P = p.parse_args()

    indir = Path(P.indir).expanduser()
    outdir = Path(P.outdir).expanduser() if P.outdir else indir

    if indir.is_file() and indir.suffix == ".dat":
        infiles = [indir]
    elif indir.is_dir():
        infiles = sorted(indir.glob("*.dat"))
    else:
        raise FileNotFoundError(f"{indir} is not a .dat file or directory")
    if not infiles:
        raise FileNotFoundError(f"no files to convert in {indir}")

    for infile in infiles:
        outfile = outdir / (infile.stem + ".h5")
        print(infile, "=>", outfile)
        dat2hdf(infile, outfile)

        if P.delete and infile != infiles[-1]:
            print("deleting", infile)
            infile.unlink()
