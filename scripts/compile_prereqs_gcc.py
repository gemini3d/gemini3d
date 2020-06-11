#!/usr/bin/env python3
"""
Installs prereqs for Gemini program for Gfortran
"""

from pathlib import Path
from argparse import ArgumentParser

import gemini3d.compile_prereqs as gc


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument(
        "libs",
        help="libraries to compile",
        choices=["netcdf", "openmpi", "hdf5", "lapack", "scalapack", "mumps"],
        nargs="+",
    )
    p.add_argument("-prefix", help="toplevel path to install libraries under", default="~/lib_gcc")
    p.add_argument("-workdir", help="toplevel path to where you keep code repos", default="~/code")
    p.add_argument("-wipe", help="wipe before completely recompiling libs", action="store_true")
    p.add_argument("-buildsys", help="build system (meson or cmake)", default="cmake")
    P = p.parse_args()

    dirs = {
        "prefix": Path(P.prefix).expanduser().resolve(),
        "workdir": Path(P.workdir).expanduser().resolve(),
    }

    # Note: HDF5 needs to be before NetCDF
    if "hdf5" in P.libs:
        gc.hdf5(dirs)
    if "netcdf" in P.libs:
        gc.netcdf_c(dirs, wipe=P.wipe)
        gc.netcdf_fortran(dirs, wipe=P.wipe)

    # Note: OpenMPI needs to be before scalapack and mumps
    if "openmpi" in P.libs:
        gc.openmpi(dirs)
    if "lapack" in P.libs:
        gc.lapack(P.wipe, dirs, P.buildsys)
    if "scalapack" in P.libs:
        gc.scalapack(P.wipe, dirs, P.buildsys)
    if "mumps" in P.libs:
        gc.mumps(P.wipe, dirs, P.buildsys)
