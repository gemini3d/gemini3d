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


def mumps(wipe: bool, dirs: typing.Dict[str, Path], buildsys: str, args: typing.List[str]):
    install_dir = dirs["prefix"] / MUMPSDIR
    source_dir = dirs["workdir"] / MUMPSDIR
    build_dir = source_dir / BUILDDIR

    update(source_dir, MUMPSGIT)

    env = get_compilers()

    if buildsys == "cmake":
        args += [f"-DCMAKE_INSTALL_PREFIX={install_dir}"]
        cmake_build(args, source_dir, build_dir, wipe, env)
    elif buildsys == "meson":
        args += [f"--prefix={dirs['prefix']}"]
        meson_build(args, source_dir, build_dir, wipe, env)
    else:
        raise ValueError(f"unknown build system {buildsys}")


@lru_cache()
def get_compilers() -> typing.Dict[str, str]:

    FC = shutil.which("ifort")
    if not FC:
        raise FileNotFoundError("Intel Fortran compiler ifort not found")
    CC = 'icl' if os.name == 'nt' else 'icc'
    CC = shutil.which(CC)
    if not CC:
        raise FileNotFoundError("Intel C compiler not found")
    CXX = 'icl' if os.name == 'nt' else 'icpc'
    CXX = shutil.which(CXX)
    if not CXX:
        raise FileNotFoundError("Intel C++ compiler not found")

    env = os.environ.update({"FC": FC, "CC": CC, "CXX": CXX})
    return env


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument(
        "-prefix", help="toplevel path to install libraries under", default="~/lib_intel"
    )
    p.add_argument(
        "-workdir", help="toplevel path to where you keep code repos", default="~/code"
    )
    p.add_argument(
        "-wipe", help="wipe before completely recompiling libs", action="store_true"
    )
    p.add_argument(
        "-b", "--buildsys", help="build system (meson or cmake)", default="cmake"
    )
    p.add_argument("-a", "--args", help="-Dfoo flags to pass to CMake / Meson", nargs="+", default=[])
    P = p.parse_args()

    dirs = {
        "prefix": Path(P.prefix).expanduser().resolve(),
        "workdir": Path(P.workdir).expanduser().resolve(),
    }

    mumps(P.wipe, dirs, P.buildsys, P.args)
