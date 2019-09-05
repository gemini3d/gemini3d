#!/usr/bin/env python3
"""
compress each file in a directory to an individual compressed file, without relative path
This allows easy transparent use of compressed files
"""
from pathlib import Path
import zipfile
import argparse

COMP = zipfile.ZIP_LZMA


def zipper(path: Path, suffix: str, overwrite: bool = False, verbose: bool = True):
    """
    Parameters
    ----------

    path: pathlib.Path
        directory to zip files under
    suffix: str
        file suffix to compress including .
    """
    path = Path(path).expanduser().resolve(True)

    flist = [f for f in path.iterdir() if f.is_file() and f.suffix == suffix]

    for file in flist:
        zipfn = file.with_suffix(".zip")
        if not overwrite and zipfn.is_file():
            try:  # rudimentary check of zip file
                with zipfile.ZipFile(zipfn, "r") as z:
                    if len(z.namelist()) > 0:
                        print("SKIP:", zipfn)
                        continue
            except zipfile.BadZipFile:
                pass

        with zipfile.ZipFile(zipfn, mode="w", compression=COMP) as z:
            print("compressing", zipfn)
            # arcname= is necessary to avoid writing vestigial (for this application) relative paths
            z.write(file, arcname=file.name)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("directory", help="directory to compress files individually")
    p.add_argument("suffix", help="file suffix to compress (including .)")
    P = p.parse_args()

    zipper(P.directory, P.suffix)
