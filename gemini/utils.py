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

from .raw import get_simsize as get_simsize_raw
from .hdf import get_simsize as get_simsize_h5

git = shutil.which("git")

Pathlike = T.Union[str, Path]

__all__ = ["gitrev", "get_cpu_count", "get_mpi_count", "get_simsize"]


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

    max_cpu = get_cpu_count(force)

    size = get_simsize(path)

    mpi_count = 1
    if size[2] == 1:  # 2D sim
        for i in range(max_cpu, 2, -1):
            mpi_count = max(math.gcd(size[1] // 2, i), mpi_count)
            if i < mpi_count:
                break
    else:  # 3D sim
        for i in range(max_cpu, 2, -1):
            mpi_count = max(math.gcd(size[2] // 2, i), mpi_count)
            if i < mpi_count:
                break

    return max(mpi_count, 1)


def get_simsize(path: Pathlike) -> T.Tuple[int, ...]:

    path = Path(path).expanduser()
    if path.is_dir():
        for suffix in (".h5", ".nc", ".dat"):
            fn = path / ("simsize" + suffix)
            if fn.is_file():
                break
    else:
        fn = path
    if not fn.is_file():
        raise FileNotFoundError(path)

    if fn.suffix == '.h5':
        return get_simsize_h5(fn)
    elif fn.suffix == '.nc':
        raise ValueError('TODO: implement NetCDF4')
    else:
        return get_simsize_raw(fn)


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
            raise SystemExit("ConnectionError: could not download {} due to {}".format(url, err))

    if filehash:
        if not file_checksum(outfile, filehash[0], filehash[1]):
            raise SystemExit("HashError: {}".format(outfile))


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
