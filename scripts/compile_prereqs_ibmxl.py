#!/usr/bin/env python3
"""
Installs prereqs for Gemini program for IBM XL compiler
"""

from pathlib import Path
from argparse import ArgumentParser

import gemini3d.compile_prereqs as gc


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument(
        "libs",
        help="libraries to compile",
        choices=["hdf5", "lapack", "scalapack", "mumps"],
        nargs="+",
    )
    p.add_argument("-prefix", help="toplevel path to install libraries under", default="~/lib_xl")
    p.add_argument("-workdir", help="toplevel path to where you keep code repos", default="~/code")
    p.add_argument("-wipe", help="wipe before completely recompiling libs", action="store_true")
    p.add_argument("-buildsys", help="build system (meson or cmake)", default="cmake")
    P = p.parse_args()

    dirs = {
        "prefix": Path(P.prefix).expanduser().resolve(),
        "workdir": Path(P.workdir).expanduser().resolve(),
    }

    env = gc.ibmxl_compilers()

    if "lapack" in P.libs:
        gc.lapack(P.wipe, dirs, env=env)
    if "scalapack" in P.libs:
        gc.scalapack(P.wipe, dirs, env=env)
    if "mumps" in P.libs:
        gc.mumps(P.wipe, dirs, env=env)
    if "hdf5" in P.libs:
        gc.hdf5(dirs, env=env)
