#!/usr/bin/env python3
"""
Installs prereqs for Gemini program for Intel compiler
"""
import subprocess
import shutil
from pathlib import Path
from argparse import ArgumentParser
import typing as T
import sys
import os
from functools import lru_cache
import gemini3d
from compile_prereqs_gcc import meson_build, cmake_build, update

# ========= user parameters ======================
BUILDDIR = "build"
NJOBS = gemini3d.get_cpu_count()

# Library parameters
HDF5VERSION = "1.10.6"
HDF5URL = "https://zenodo.org/record/3659270/files/hdf5-1.10.6.tar.bz2?download=1"
HDF5MD5 = "03095102a6118c32a75a9b9b40be66f2"
HDF5DIR = f"hdf5-{HDF5VERSION}"

MUMPSGIT = "https://github.com/scivision/mumps"
MUMPSDIR = "mumps"

# ========= end of user parameters ================

nice = ["nice"] if sys.platform == "linux" else []


def hdf5(dirs: T.Dict[str, Path]):
    """ build and install HDF5 """
    if os.name == "nt":
        raise SystemExit("Please use binaries from HDF Group for Windows appropriate for your compiler.")

    install_dir = dirs["prefix"] / HDF5DIR
    source_dir = dirs["workdir"] / HDF5DIR

    tarfn = f"hdf5-{HDF5VERSION}.tar.bz2"
    gemini3d.url_retrieve(HDF5URL, tarfn, ("md5", HDF5MD5))
    gemini3d.extract_tar(tarfn, source_dir)

    env = get_compilers()

    cmd = nice + ["./configure", f"--prefix={install_dir}", "--enable-fortran", "--enable-build-mode=production"]

    subprocess.check_call(cmd, cwd=source_dir, env=env)

    cmd = nice + ["make", "-C", str(source_dir), f"-j {NJOBS}", "install"]
    subprocess.check_call(cmd)


def mumps(wipe: bool, dirs: T.Dict[str, Path], buildsys: str):
    """ build and install Mumps """
    install_dir = dirs["prefix"] / MUMPSDIR
    source_dir = dirs["workdir"] / MUMPSDIR
    build_dir = source_dir / BUILDDIR

    update(source_dir, MUMPSGIT)

    if buildsys == "cmake":
        args = [f"-DCMAKE_INSTALL_PREFIX={install_dir}"]
        cmake_build(args, source_dir, build_dir, wipe, env=get_compilers())
    elif buildsys == "meson":
        args = [f"--prefix={dirs['prefix']}"]
        meson_build(args, source_dir, build_dir, wipe, env=get_compilers())
    else:
        raise ValueError(f"unknown build system {buildsys}")


@lru_cache()
def get_compilers() -> T.Mapping[str, str]:
    """ get paths to Intel compilers """
    env = os.environ

    fc_name = "ifort"
    cc_name = "icl" if os.name == "nt" else "icc"
    cxx_name = "icl" if os.name == "nt" else "icpc"

    fc = env.get("FC", "")
    if fc_name not in fc:
        fc = shutil.which(fc_name)
    if not fc:
        raise FileNotFoundError(f"{fc_name} not found")

    cc = env.get("CC", "")
    if cc_name not in cc:
        cc = shutil.which(cc_name)
    if not cc:
        raise FileNotFoundError(f"{cc_name} not found")

    cxx = env.get("CXX", "")
    if cxx_name not in cxx:
        cxx = shutil.which(cxx_name)
    if not cxx:
        raise FileNotFoundError(f"{cxx_name} not found")

    env.update({"FC": fc, "CC": cc, "CXX": cxx})

    return env


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("libs", help="libraries to compile", choices=["hdf5", "mumps"], nargs="+")
    p.add_argument("-prefix", help="toplevel path to install libraries under", default="~/lib_intel")
    p.add_argument("-workdir", help="toplevel path to where you keep code repos", default="~/code")
    p.add_argument("-wipe", help="wipe before completely recompiling libs", action="store_true")
    p.add_argument("-b", "--buildsys", help="build system (meson or cmake)", default="cmake")
    P = p.parse_args()

    dirs = {"prefix": Path(P.prefix).expanduser().resolve(), "workdir": Path(P.workdir).expanduser().resolve()}

    if "mumps" in P.libs:
        mumps(P.wipe, dirs, P.buildsys)
    if "hdf5" in P.libs:
        hdf5(dirs)
