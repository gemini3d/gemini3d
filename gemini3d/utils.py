import subprocess
import os
import shutil
from pathlib import Path
import math
import typing as T
import hashlib
import urllib.request
import urllib.error
import socket
import zipfile
import tarfile

try:
    import psutil
except ImportError:
    psutil = None
    # pip install psutil will improve CPU utilization.

from .config import read_config
from .base import get_simsize

git = shutil.which("git")

Pathlike = T.Union[str, Path]

__all__ = ["gitrev", "get_cpu_count", "get_mpi_count"]


def gitrev() -> str:
    if not git:
        return ""

    return subprocess.check_output([git, "rev-parse", "--short", "HEAD"], universal_newlines=True).strip()


def get_cpu_count(force: int = None) -> int:
    if force:
        max_cpu = force
        extradiv = 1
    else:
        max_cpu = None
        # without psutil, hyperthreaded CPU may overestimate physical count by factor of 2 (or more)
        if psutil is not None:
            max_cpu = psutil.cpu_count(logical=False)
            extradiv = 1
            if max_cpu is None:
                max_cpu = psutil.cpu_count()
                extradiv = 2
        if max_cpu is None:
            max_cpu = os.cpu_count()
            extradiv = 2

    return max_cpu // extradiv


def get_mpi_count(path: Pathlike, force: int = None) -> int:

    path = Path(path).expanduser()

    max_cpu = get_cpu_count(force)

    # %% config.nml file or directory or simsize.h5?
    if path.is_dir():
        size = get_simsize(path)
    elif path.is_file():
        if path.suffix in (".h5", ".nc", ".dat"):
            size = get_simsize(path)
        elif path.suffix in (".ini", ".nml"):
            params = read_config(path)
            # OK to use indat_size because we're going to run a sim on this machine
            size = get_simsize(params["indat_size"])
    else:
        raise FileNotFoundError(f"{path} is not a file or directory")

    mpi_count = 1
    if size[2] == 1:
        # 2D sim
        for i in range(max_cpu, 2, -1):
            mpi_count = max(math.gcd(size[1], i), mpi_count)
            if i < mpi_count:
                break
    else:
        # 3D sim
        for i in range(max_cpu, 2, -1):
            mpi_count = max(math.gcd(size[2], i), mpi_count)
            if i < mpi_count:
                break

    return mpi_count


def url_retrieve(url: str, outfile: Pathlike, filehash: T.Sequence[str] = None, overwrite: bool = False):
    """
    Parameters
    ----------
    url: str
        URL to download from
    outfile: pathlib.Path
        output filepath (including name)
    filehash: tuple of str, str
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
            raise ConnectionError(f"could not download {url} due to {err}")

    if filehash:
        if not file_checksum(outfile, filehash[0], filehash[1]):
            raise ValueError(f"Hash mismatch: {outfile}")


def file_checksum(fn: Path, mode: str, filehash: str) -> bool:
    h = hashlib.new(mode)
    h.update(fn.read_bytes())
    return h.hexdigest() == filehash


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
        # tarfile gives confusing error on missing file
        raise FileNotFoundError(fn)
    with tarfile.open(fn) as z:
        z.extractall(str(outpath.parent))
