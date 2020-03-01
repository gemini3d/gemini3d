#!/usr/bin/env python
"""
runs a job
"""
import argparse
import gemini.job
import sys
from time import monotonic


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("config_file", help="path to config*.nml file")
    p.add_argument("out_dir", help="simulation output directory")
    p.add_argument("-mpiexec", help="path to desired mpiexec executable")
    p.add_argument("-gemexe", help="path to desired gemini.bin")
    p.add_argument("-matlab", help="Use Matlab instead of Python", action="store_true")
    p.add_argument("-f", "--force", help="force regeneration of simulation", action="store_true")
    P = p.parse_args()

    ret = -1

    params = {
        "config_file": P.config_file,
        "out_dir": P.out_dir,
        "mpiexec": P.mpiexec,
        "gemexe": P.gemexe,
        "matlab": P.matlab,
        "force": P.force,
    }

    tic = monotonic()
    try:
        ret = gemini.job.runner(params)
    except FileNotFoundError as excp:
        print(excp, file=sys.stderr)
        print(
            "\nA necessary simulation input file was not found."
            "\nThis can mean that the simulation initialization script wasn't run first.\n",
            file=sys.stderr,
        )
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
