#!/usr/bin/env python3
"""
run test
"""
import argparse
import sys
import typing
import subprocess
from pathlib import Path
import gemini

R = Path(__file__).resolve().parents[1]
Rt = R / "tests/data"

zenodo: typing.Dict[str, typing.Any] = {
    "2d_fang": {
        "url": "https://zenodo.org/record/3677638/files/test2d_fang.zip?download=1",
        "md5": "03c183bbc91706223313e5c15771918e",
        "dir": Rt / "test2d_fang",
        "zip": Rt / "test2d_fang.zip",
    },
    "2d_glow": {
        "url": "https://zenodo.org/record/3677638/files/test2d_glow.zip?download=1",
        "md5": "bd9a9c38bb462cc22cc0ea232e03dc21",
        "dir": Rt / "test2d_glow",
        "zip": Rt / "test2d_glow.zip",
    },
    "3d_fang": {
        "url": "https://zenodo.org/record/3687202/files/test3d_fang.zip?download=1",
        "md5": "b4c5fc43243b33549b8324c9a56ee198",
        "dir": Rt / "test3d_fang",
        "zip": Rt / "test3d_fang.zip",
    },
    "3d_glow": {
        "url": "https://zenodo.org/record/3687202/files/test3d_glow.zip?download=1",
        "md5": "d70c8ee5a699ae028b0ffecb750fb5c6",
        "dir": Rt / "test3d_glow",
        "zip": Rt / "test3d_glow.zip",
    },
}


def run_test(testname: str, mpiexec: str, exe: str, nml: str, outdir: str, mpi_count: int = None):
    z = zenodo[testname]
    try:
        gemini.url_retrieve(z["url"], z["zip"], ("md5", z["md5"]))
    except (ConnectionError, ValueError) as e:
        print(f"SKIP: problem downloading reference data {e}", file=sys.stderr)
        raise SystemExit(77)

    gemini.extract_zip(z["zip"], z["dir"])

    if not Path(nml).expanduser().is_file():
        print(f"SKIP: {nml} does not exist", file=sys.stderr)
        raise SystemExit(77)

    if not mpi_count:
        mpi_count = gemini.get_mpi_count(z["dir"] / "inputs/simsize.h5")

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
