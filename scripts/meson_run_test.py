#!/usr/bin/env python3
"""
run test
"""
import argparse
import sys
import typing
import subprocess
from pathlib import Path
import gemini3d.utils as gu

R = Path(__file__).resolve().parents[1]
Rt = R / "tests/data"

zenodo: typing.Dict[str, typing.Any] = {
    "2d_fang": {
        "url": "https://zenodo.org/record/3697235/files/test2d_fang.zip?download=1",
        "md5": "7b237883b17d6da331b7db0309cc3da4",
        "dir": Rt / "test2d_fang",
        "zip": Rt / "test2d_fang.zip",
    },
    "2d_glow": {
        "url": "https://zenodo.org/record/3697235/files/test2d_glow.zip?download=1",
        "md5": "4b8e31e239918ff9e02cbc3033099d60",
        "dir": Rt / "test2d_glow",
        "zip": Rt / "test2d_glow.zip",
    },
    "3d_fang": {
        "url": "https://zenodo.org/record/3698137/files/test3d_fang.zip?download=1",
        "md5": "a270586fbd41ff58e3a5ed0acac63708",
        "dir": Rt / "test3d_fang",
        "zip": Rt / "test3d_fang.zip",
    },
    "3d_glow": {
        "url": "https://zenodo.org/record/3698137/files/test3d_glow.zip?download=1",
        "md5": "993d817ee847314a631396608df56e0c",
        "dir": Rt / "test3d_glow",
        "zip": Rt / "test3d_glow.zip",
    },
}


def run_test(testname: str, mpiexec: str, exe: str, nml: str, outdir: str, mpi_count: int = None):
    z = zenodo[testname]
    try:
        gu.url_retrieve(z["url"], z["zip"], ("md5", z["md5"]))
    except (ConnectionError, ValueError) as e:
        print(f"SKIP: problem downloading reference data {e}", file=sys.stderr)
        raise SystemExit(77)

    gu.extract_zip(z["zip"], z["dir"])

    if not Path(nml).expanduser().is_file():
        print(f"SKIP: {nml} does not exist", file=sys.stderr)
        raise SystemExit(77)

    if not mpi_count:
        mpi_count = gu.get_mpi_count(z["dir"] / "inputs/simsize.h5")

    # have to get exe as absolute path
    exe_abs = Path(exe).resolve()
    cmd = [mpiexec, "-np", str(mpi_count), str(exe_abs), nml, outdir]
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
    P = p.parse_args()

    run_test(P.testname, P.mpiexec, P.exe, P.nml, P.outdir, mpi_count=P.np)
