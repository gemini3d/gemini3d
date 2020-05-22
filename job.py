#!/usr/bin/env python
"""
runs a job
"""
import argparse
from pathlib import Path
import gemini3d.job
import sys
from time import monotonic


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("config_file", help="path to config*.nml file")
    p.add_argument("out_dir", help="simulation output directory")
    p.add_argument("-mpiexec", help="path to desired mpiexec executable")
    p.add_argument("-gemexe", help="path to desired gemini.bin")
    p.add_argument("-n", "--cpu", help="number of CPU cores", type=int, default=0)
    p.add_argument("-f", "--force", help="force regeneration of simulation", action="store_true")
    p.add_argument("-out_format", help="override Fortran output file format", choices=["h5", "nc", "raw"])
    P = p.parse_args()

    ret = -1

    params = {
        "config_file": Path(P.config_file).expanduser(),
        "out_dir": Path(P.out_dir).expanduser().resolve(),
        "mpiexec": P.mpiexec,
        "gemexe": P.gemexe,
        "force": P.force,
        "out_format": P.out_format,
        "cpu_count": P.cpu,
    }

    tic = monotonic()
    try:
        ret = gemini3d.job.runner(params)
    except FileNotFoundError:
        print(
            "\nA necessary simulation input file was not found."
            "\nThis can mean that the simulation initialization script wasn't run first.\n",
            file=sys.stderr,
        )
        raise
    except EnvironmentError as excp:
        print(excp, file=sys.stderr)
        print(
            "If you need to build Gemini, from the Gemini directory:\n\n",
            "cmake -B build\n",
            "cmake --build build\n",
            file=sys.stderr,
        )

    print(f"job.py ran in {monotonic() - tic:.3f} seconds.")

    raise SystemExit(ret)
