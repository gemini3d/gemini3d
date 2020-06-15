#!/usr/bin/env python3
"""
Installs prereqs for Gemini program for Intel compiler
"""

from pathlib import Path
from argparse import ArgumentParser

import gemini3d.compile_prereqs as gc


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument(
        "libs", help="libraries to compile", choices=["netcdf", "hdf5", "mumps"], nargs="+"
    )
    p.add_argument(
        "-prefix", help="toplevel path to install libraries under", default="~/lib_intel"
    )
    p.add_argument("-workdir", help="toplevel path to where you keep code repos", default="~/code")
    p.add_argument("-wipe", help="wipe before completely recompiling libs", action="store_true")
    P = p.parse_args()

    dirs = {
        "prefix": Path(P.prefix).expanduser().resolve(),
        "workdir": Path(P.workdir).expanduser().resolve(),
    }

    env = gc.intel_compilers()
    if "hdf5" in P.libs:
        gc.hdf5(dirs, env=env)
    if "netcdf" in P.libs:
        gc.netcdf_c(dirs, env=env, wipe=P.wipe)
        gc.netcdf_fortran(dirs, env=env, wipe=P.wipe)

    if "mumps" in P.libs:
        gc.mumps(P.wipe, dirs, env=env)
