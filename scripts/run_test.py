#!/usr/bin/env python3
"""
run test
"""
import argparse
from configparser import ConfigParser
import sys
import typing
import subprocess
from pathlib import Path
import typing as T
import logging

from web import url_retrieve, extract_zip

try:
    from gemini3d import get_mpi_count
except ImportError as e:
    logging.warning(f"could not use Gemini3D get_mpi_count(), falling back to single CPU.   {e}")

    def get_mpi_count(path: Path) -> int:
        return 1


R = Path(__file__).resolve().parents[1]
Rt = R / "tests/data"

Pathlike = T.Union[Path, str]


def get_test_params(test_name: str, url_file: Path) -> T.Dict[str, T.Any]:
    """ get URL and MD5 for a test name """
    ini = Path(url_file).expanduser().read_text()
    C = ConfigParser()
    C.read_string(ini)

    z = {
        "url": C.get(test_name, "url"),
        "md5": C.get(test_name, "md5"),
        "dir": Rt / f"test{test_name}",
        "zip": Rt / f"test{test_name}.zip",
    }

    return z


def download_and_extract(z: T.Dict[str, T.Any], url_ini: Path):

    try:
        url_retrieve(z["url"], z["zip"], ("md5", z["md5"]))
    except (ConnectionError, ValueError) as e:
        print(f"SKIP: problem downloading reference data {e}", file=sys.stderr)
        raise SystemExit(77)
    except KeyError as e:
        print(f"SKIP: problem getting reference config from {url_ini} {e}", file=sys.stderr)
        raise SystemExit(77)


def run_test(testname: str, mpiexec: str, exe: str, nml: str, outdir: str, mpi_count: int = None, out_format: str = None):
    """ configure and run a test
    This is usually called from CMake Ctest
    """

    url_ini = R / "gemini3d/tests/url.ini"

    z = get_test_params(testname, url_ini)

    if not z["dir"].is_dir():
        download_and_extract(z, url_ini)

    extract_zip(z["zip"], z["dir"])

    if not Path(nml).expanduser().is_file():
        print(f"SKIP: file {nml} does not exist", file=sys.stderr)
        raise SystemExit(77)

    if not mpi_count:
        mpi_count = get_mpi_count(z["dir"] / "inputs/simsize.h5")

    # have to get exe as absolute path
    exe_abs = Path(exe).resolve()
    cmd = [mpiexec, "-np", str(mpi_count), str(exe_abs), nml, outdir]
    if out_format:
        cmd += ["-out_format", out_format]
    print(" ".join(cmd))
    ret = subprocess.run(cmd, cwd=R)

    raise SystemExit(ret.returncode)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("testname")
    p.add_argument("mpiexec")
    p.add_argument("exe")
    p.add_argument("nml")
    p.add_argument("outdir")
    p.add_argument("-np", help="force number of MPI images", type=int)
    p.add_argument("-out_format", help="override config.nml output file format", choices=["h5", "nc", "raw"])
    P = p.parse_args()

    run_test(P.testname, P.mpiexec, P.exe, P.nml, P.outdir, mpi_count=P.np, out_format=P.out_format)
