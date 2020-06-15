#!/usr/bin/env python
"""
CLI scripts
"""
import argparse
from pathlib import Path
import sys
from time import monotonic

from .job import runner
from .compare import compare_all
from .model_setup import model_setup


def gemini_setup():
    p = argparse.ArgumentParser()
    p.add_argument("config_file", help="path to config*.nml file")
    p.add_argument("out_dir", help="simulation output directory")
    P = p.parse_args()

    model_setup(P.config_file, P.out_dir)


def compare_run():
    tol = {
        "rtol": 1e-5,
        "rtolN": 1e-5,
        "rtolT": 1e-5,
        "rtolJ": 1e-5,
        "rtolV": 1e-5,
        "atol": 1e-8,
        "atolN": 1e9,
        "atolT": 100,
        "atolJ": 1e-7,
        "atolV": 50,
    }

    p = argparse.ArgumentParser(description="Compare simulation file outputs and inputs")
    p.add_argument("outdir", help="directory to compare")
    p.add_argument("refdir", help="reference directory")
    p.add_argument("-p", "--plot", help="make plots of differences", action="store_true")
    p.add_argument("-only", help="only check in or out", choices=["in", "out"])
    p.add_argument(
        "-file_format",
        help="specify file format to read from output dir",
        choices=["h5", "nc", "raw"],
    )
    P = p.parse_args()

    out_errs, in_errs = compare_all(P.outdir, P.refdir, tol, P.plot, P.file_format, P.only)

    if out_errs or in_errs:
        raise SystemExit(f"{out_errs} output errors, {in_errs} input errors")
    else:
        only = ("out", "in") if not P.only else P.only

        print(f"OK: Gemini {only} comparison {P.outdir} {P.refdir}")


def gemini_run():
    p = argparse.ArgumentParser()
    p.add_argument("config_file", help="path to config*.nml file")
    p.add_argument("out_dir", help="simulation output directory")
    p.add_argument("-mpiexec", help="path to desired mpiexec executable")
    p.add_argument("-gemexe", help="path to desired gemini.bin")
    p.add_argument("-n", "--cpu", help="number of CPU cores", type=int, default=0)
    p.add_argument("-f", "--force", help="force regeneration of simulation", action="store_true")
    p.add_argument(
        "-out_format", help="override Fortran output file format", choices=["h5", "nc", "raw"]
    )
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
        ret = runner(params)
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
