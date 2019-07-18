#!/usr/bin/env python
from pathlib import Path
import argparse
import zipfile


def extract_files(zipfn: Path, outpath: Path, overwrite: bool = False):
    outpath = Path(outpath).expanduser().resolve()
    # need .resolve() in case intermediate relative dir doesn't exist
    if outpath.is_dir() and not overwrite:
        return
    zipfn = Path(zipfn).expanduser().resolve()
    with zipfile.ZipFile(zipfn) as z:
        z.extractall(str(outpath))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('zipfile', help='.zip file to extract')
    p.add_argument('outpath', help='path to extract into')
    P = p.parse_args()

    extract_files(P.zipfile, P.outpath)
