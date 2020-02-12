#!/usr/bin/env python3
from pathlib import Path
import argparse
import zipfile
import tarfile
import typing as T

Pathlike = T.Union[str, Path]


def extract_zip(fn: Pathlike, outpath: Pathlike, overwrite: bool = False):
    outpath = Path(outpath).expanduser().resolve()
    # need .resolve() in case intermediate relative dir doesn't exist
    if outpath.is_dir() and not overwrite:
        return

    fn = Path(fn).expanduser().resolve()
    with zipfile.ZipFile(fn) as z:
        z.extractall(str(outpath.parent))


def extract_tar(fn: Pathlike, outpath: Pathlike, overwrite: bool = False):
    outpath = Path(outpath).expanduser().resolve()
    # need .resolve() in case intermediate relative dir doesn't exist
    if outpath.is_dir() and not overwrite:
        return

    fn = Path(fn).expanduser().resolve()
    if not fn.is_file():
        raise FileNotFoundError(fn)  # keep this, tarfile gives confusing error
    with tarfile.open(fn) as z:
        z.extractall(str(outpath.parent))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("infile", help="compressed file to extract")
    p.add_argument("outpath", help="path to extract into")
    P = p.parse_args()

    infile = Path(P.infile)
    if infile.suffix.lower() == ".zip":
        extract_zip(infile, P.outpath)
    elif infile.suffix.lower() in (".tar", ".gz", ".bz2", ".xz"):
        extract_tar(infile, P.outpath)
    else:
        raise ValueError("Not sure how to decompress {}".format(infile))
