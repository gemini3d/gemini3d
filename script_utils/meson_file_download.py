#!/usr/bin/env python3
"""
We use SystemExit as this will not blast the whole traceback to Meson.
Usually just a terse stderr will suffice and not overwhelm the Meson user.
"""
from pathlib import Path
import urllib.request
import urllib.error
import hashlib
import argparse
import typing
import socket


def url_retrieve(
    url: str, outfile: Path, hash: typing.Sequence[str] = None, overwrite: bool = False
):
    """
    Parameters
    ----------
    url: str
        URL to download from
    outfile: pathlib.Path
        output filepath (including name)
    hash: tuple of str, str
        hash type (md5, sha1, etc.) and hash
    overwrite: bool
        overwrite if file exists
    """
    outfile = Path(outfile).expanduser().resolve()
    if outfile.is_dir():
        raise ValueError("Please specify full filepath, including filename")
    # need .resolve() in case intermediate relative dir doesn't exist
    if overwrite or not outfile.is_file():
        outfile.parent.mkdir(parents=True, exist_ok=True)
        try:
            urllib.request.urlretrieve(url, str(outfile))
        except (socket.gaierror, urllib.error.URLError) as err:
            raise SystemExit(
                "ConnectionError: could not download {} due to {}".format(url, err)
            )

    if hash:
        if not file_checksum(outfile, hash[0], hash[1]):
            raise SystemExit("HashError: {}".format(outfile))


def file_checksum(fn: Path, mode: str, hash: str) -> bool:
    h = hashlib.new(mode)
    h.update(fn.read_bytes())
    return h.hexdigest() == hash


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("url", help="URL to file download")
    p.add_argument("outfile", help="filename to download to")
    p.add_argument("-hash", help="expected hash", nargs=2)
    P = p.parse_args()

    url_retrieve(P.url, P.outfile, P.hash)
