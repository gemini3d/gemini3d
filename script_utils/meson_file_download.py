#!/usr/bin/env python
from pathlib import Path
import requests
import hashlib
import argparse


def url_retrieve(url: str, outfile: Path, md5sum: str = None, overwrite: bool = False):
    outfile = Path(outfile).expanduser().resolve()
    # need .resolve() in case intermediate relative dir doesn't exist
    if outfile.is_file() and not overwrite:
        return
    outfile.parent.mkdir(parents=True, exist_ok=True)

    R = requests.get(url, allow_redirects=True)
    if R.status_code != 200:
        raise ConnectionError('could not download {}\nerror code: {}'.format(url, R.status_code))

    outfile.write_bytes(R.content)

    if md5sum:
        if not file_checksum(outfile, md5sum, 'md5'):
            raise OSError('hash mismatch: Failed to download {}'.format(outfile))


def file_checksum(fn: Path, hash: str, mode: str) -> bool:
    h = hashlib.new(mode)
    h.update(fn.read_bytes())
    return h.hexdigest() == hash


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('url', help='URL to file download')
    p.add_argument('outfile', help='filename to download to')
    p.add_argument('md5sum', help='expected MD5 hash', nargs='?')
    P = p.parse_args()

    url_retrieve(P.url, P.outfile, P.md5sum)
