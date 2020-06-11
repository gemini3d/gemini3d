#!/usr/bin/env python3
"""
Installs prereqs for Gemini program for Intel compiler
"""

import shutil
from pathlib import Path
from argparse import ArgumentParser
import typing as T
import sys
import os

import gemini3d.compile_prereqs as gc


nice = ["nice"] if sys.platform == "linux" else []


def get_compilers() -> T.Mapping[str, str]:
    """ get paths to compilers """
    env = os.environ

    fc_name = "ifort"
    cc_name = "icl" if os.name == "nt" else "icc"
    cxx_name = "icl" if os.name == "nt" else "icpc"

    fc = env.get("FC", "")
    if fc_name not in fc:
        fc = shutil.which(fc_name)
    if not fc:
        raise FileNotFoundError(fc_name)

    cc = env.get("CC", "")
    if cc_name not in cc:
        cc = shutil.which(cc_name)
    if not cc:
        raise FileNotFoundError(cc_name)

    cxx = env.get("CXX", "")
    if cxx_name not in cxx:
        cxx = shutil.which(cxx_name)
    if not cxx:
        raise FileNotFoundError(cxx_name)

    env.update({"FC": fc, "CC": cc, "CXX": cxx})

    return env


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
    p.add_argument("-buildsys", help="build system (meson or cmake)", default="cmake")
    P = p.parse_args()

    dirs = {
        "prefix": Path(P.prefix).expanduser().resolve(),
        "workdir": Path(P.workdir).expanduser().resolve(),
    }

    env = get_compilers()
    if "hdf5" in P.libs:
        gc.hdf5(dirs, env=env)
    if "netcdf" in P.libs:
        gc.netcdf_c(dirs, env=env, wipe=P.wipe)
        gc.netcdf_fortran(dirs, env=env, wipe=P.wipe)

    if "mumps" in P.libs:
        gc.mumps(P.wipe, dirs, P.buildsys, env=env)
