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
import shutil
import importlib.resources

import gemini3d

Rtop = Path(__file__).resolve().parents[1]
Rtest = Rtop / "tests/data"


def get_test_params(test_name: str, url_file: Path) -> T.Dict[str, T.Any]:
    """ get URL and MD5 for a test name """
    ini = Path(url_file).expanduser().read_text()
    C = ConfigParser()
    C.read_string(ini)

    z = {
        "url": C.get(test_name, "url"),
        "md5": C.get(test_name, "md5"),
        "dir": Rtest / f"test{test_name}",
        "zip": Rtest / f"test{test_name}.zip",
    }

    return z


def download_and_extract(z: T.Dict[str, T.Any], url_ini: Path):

    try:
        gemini3d.url_retrieve(z["url"], z["zip"], ("md5", z["md5"]))
    except (ConnectionError, ValueError) as e:
        print(f"SKIP: problem downloading reference data {e}", file=sys.stderr)
        raise SystemExit(77)
    except KeyError as e:
        print(
            f"SKIP: problem getting reference config from {url_ini} {e}", file=sys.stderr,
        )
        raise SystemExit(77)


def run_test(
    testname: str,
    mpiexec: str,
    exe: str,
    outdir: Path,
    *,
    mpi_count: int = None,
    out_format: str = None,
    dryrun: bool = False,
) -> int:
    """ configure and run a test
    This is usually called from CMake Ctest
    """

    with importlib.resources.path("gemini3d.tests", "gemini3d_url.ini") as url_ini:
        z = get_test_params(testname, url_ini)

        if not z["dir"].is_dir():
            download_and_extract(z, url_ini)

    gemini3d.extract_zip(z["zip"], z["dir"])

    outdir = Path(outdir).expanduser().resolve()
    # prepare simulation output directory
    input_dir = outdir / "inputs"
    nml = z["dir"] / "inputs/config.nml"
    input_dir.mkdir(parents=True, exist_ok=True)

    # a normal, non-test simulation already has all these files in the
    # output directory. here, we use a reference simulation input data.
    # the input data generation is tested elsewhere in PyGemini.
    # Here, we want to test that we can create
    # data match to reference outputs from reference inputs.
    shutil.copy2(nml, input_dir)
    cfg = gemini3d.read_config(nml)
    shutil.copy2(z["dir"] / cfg["indat_size"], input_dir)
    shutil.copy2(z["dir"] / cfg["indat_grid"], input_dir)
    shutil.copy2(z["dir"] / cfg["indat_file"], input_dir)
    if "precdir" in cfg and not (outdir / cfg["precdir"]).is_dir():
        shutil.copytree(z["dir"] / cfg["precdir"], outdir / cfg["precdir"])
    if "E0dir" in cfg and not (outdir / cfg["E0dir"]).is_dir():
        shutil.copytree(z["dir"] / cfg["E0dir"], outdir / cfg["E0dir"])
    if "neutral_perturb" in cfg and not (outdir / cfg["sourcedir"]).is_dir():
        shutil.copytree(z["dir"] / cfg["sourcedir"], outdir / cfg["sourcedir"])

    if not mpi_count:
        mpi_count = gemini3d.get_mpi_count(z["dir"] / cfg["indat_size"], 0)

    # have to get exe as absolute path
    exe_abs = Path(exe).resolve()
    if mpiexec:
        cmd = [mpiexec, "-np", str(mpi_count), str(exe_abs), str(outdir)]
    else:
        cmd = [str(exe_abs), str(outdir)]
    if out_format:
        cmd += ["-out_format", out_format]
    if dryrun:
        cmd.append("-dryrun")
    print(" ".join(cmd))
    ret = subprocess.run(cmd, cwd=Rtop)
    return ret.returncode


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("testname")
    p.add_argument("-mpiexec", help="mpiexec path")
    p.add_argument("exe")
    p.add_argument("outdir")
    p.add_argument("-np", help="force number of MPI images", type=int)
    p.add_argument(
        "-out_format", help="override config.nml output file format", choices=["h5", "nc", "raw"],
    )
    p.add_argument("-dryrun", help="only run first time step", action="store_true")
    P = p.parse_args()

    ret = run_test(
        P.testname,
        P.mpiexec,
        P.exe,
        P.outdir,
        mpi_count=P.np,
        out_format=P.out_format,
        dryrun=P.dryrun,
    )

    raise SystemExit(ret)
