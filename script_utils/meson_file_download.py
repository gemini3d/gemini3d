#!/usr/bin/env python
from pathlib import Path
import urllib.request
import hashlib
import argparse
import typing


def url_retrieve(url: str,
                 outfile: Path,
                 hash: typing.Sequence[str] = None,
                 overwrite: bool = False):
    """
    Parameters
    ----------
    url: str
        URL to downlaod from
    outfile: pathlib.Path
        output filepath (including name)
    hash: tuple of str, str
        hash type (md5, sha1, etc.) and hash
    overwrite: bool
        overwrite if file exists
    """
    outfile = Path(outfile).expanduser().resolve()
    # need .resolve() in case intermediate relative dir doesn't exist
    if outfile.is_file() and not outfile.is_dir() and not overwrite:
        return
    outfile.parent.mkdir(parents=True, exist_ok=True)

    urllib.request.urlretrieve(url, str(outfile))

    if hash:
        if not file_checksum(outfile, hash[0], hash[1]):
            raise OSError('hash mismatch: Failed to download {}'.format(outfile))


def file_checksum(fn: Path, mode: str, hash: str) -> bool:
    h = hashlib.new(mode)
    h.update(fn.read_bytes())
    return h.hexdigest() == hash


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('url', help='URL to file download')
    p.add_argument('outfile', help='filename to download to')
    p.add_argument('-hash', help='expected hash', nargs='2')
    P = p.parse_args()

    url_retrieve(P.url, P.outfile, P.hash)
