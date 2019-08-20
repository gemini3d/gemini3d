#!/usr/bin/env python3
"""
Installs prereqs for Gemini program for Intel compiler
"""
import shutil
from pathlib import Path
from argparse import ArgumentParser
import typing
import sys
import os
from functools import lru_cache

from compile_prereqs_gcc import meson_build, cmake_build, update

# ========= user parameters ======================
BUILDDIR = "build"
LOADLIMIT = "-l 4"

# Library parameters

MUMPSGIT = "https://github.com/scivision/mumps"
MUMPSDIR = "mumps"

# ========= end of user parameters ================

nice = ["nice"] if sys.platform == "linux" else []


def mumps(wipe: bool, dirs: typing.Dict[str, Path], buildsys: str):
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
def get_compilers() -> typing.Mapping[str, str]:

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
    p.add_argument(
        "-prefix", help="toplevel path to install libraries under", default="~/lib_intel"
    )
    p.add_argument("-workdir", help="toplevel path to where you keep code repos", default="~/code")
    p.add_argument("-wipe", help="wipe before completely recompiling libs", action="store_true")
    p.add_argument("-b", "--buildsys", help="build system (meson or cmake)", default="cmake")
    P = p.parse_args()

    dirs = {
        "prefix": Path(P.prefix).expanduser().resolve(),
        "workdir": Path(P.workdir).expanduser().resolve(),
    }

    mumps(P.wipe, dirs, P.buildsys)
