#!/usr/bin/env python
"""
Installs prereqs for Gemini program for Intel compiler
"""
import subprocess
import shutil
import logging
from pathlib import Path
import pkg_resources
from argparse import ArgumentParser
import typing
import sys
import os
from functools import lru_cache

# ========= user parameters ======================
BUILDDIR = "build"
LOADLIMIT = "-l 4"

# Library parameters

MUMPSGIT = "https://github.com/scivision/mumps"
MUMPSDIR = "mumps"

# ========= end of user parameters ================

nice = ["nice"] if sys.platform == "linux" else []

ENV = os.environ


def mumps(wipe: bool, dirs: typing.Dict[str, Path], buildsys: str):
    install_lib = dirs["prefix"] / MUMPSDIR
    source_lib = dirs["workdir"] / MUMPSDIR
    build_lib = source_lib / BUILDDIR

    update(source_lib, MUMPSGIT)

    if buildsys == "cmake":
        args = [
            f"-DCMAKE_INSTALL_PREFIX={install_lib}" "-B",
            str(build_lib),
            "-S",
            str(source_lib),
        ]
        cmake_build(args, build_lib, wipe)
    elif buildsys == "meson":
        args = [f"--prefix={dirs['prefix']}"]
        meson_build(args, build_lib, wipe)
    else:
        raise ValueError(f"unknown build system {buildsys}")


def cmake_build(args: typing.List[str], build_dir: Path, wipe: bool):
    cmake = cmake_minimum_version("3.13")
    cachefile = build_dir / "CMakeCache.txt"
    if wipe and cachefile.is_file():
        cachefile.unlink()

    subprocess.check_call(nice + [cmake] + args, env=ENV)

    subprocess.check_call(
        nice + [cmake, "--build", str(build_dir), "--parallel", "--target", "install"]
    )

    subprocess.check_call(
        nice + ["ctest", "--parallel", "--output-on-failure"], cwd=str(build_dir)
    )


def meson_build(args: typing.List[str], build_dir: Path, wipe: bool):
    meson = shutil.which("meson")
    if not meson:
        raise FileNotFoundError("Meson not found.")

    if wipe and (build_dir / "ninja.build").is_file():
        args.append("--wipe")

    subprocess.check_call(nice + [meson, "setup"] + args + [str(build_dir)], env=ENV)

    for op in ("test", "install"):
        subprocess.check_call(nice + [meson, op, "-C", str(build_dir)])


@lru_cache()
def cmake_minimum_version(min_version: str = None) -> str:
    """
    if CMake is at least minimum version, return path to CMake executable
    """

    cmake = shutil.which("cmake")
    if not cmake:
        raise FileNotFoundError("could not find CMake")

    if not min_version:
        return cmake
    cmake_ver = (
        subprocess.check_output([cmake, "--version"], universal_newlines=True)
        .split("\n")[0]
        .split(" ")[2]
    )
    if pkg_resources.parse_version(cmake_ver) < pkg_resources.parse_version(
        min_version
    ):
        raise ValueError(
            f"CMake {cmake_ver} is less than minimum required {min_version}"
        )

    return cmake


def update(path: Path, repo: str):
    """
    Use Git to update a local repo, or clone it if not already existing.

    we use cwd= instead of "git -C" for very old Git versions that might be on your HPC.
    """
    GITEXE = shutil.which("git")

    if not GITEXE:
        logging.warning("Git not available, cannot check for library updates")
        return

    if path.is_dir():
        subprocess.check_call([GITEXE, "pull"], cwd=str(path))
    else:
        subprocess.check_call([GITEXE, "clone", repo, str(path)])


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument(
        "-prefix", help="toplevel path to install libraries under", default="~/lib_gcc"
    )
    p.add_argument(
        "-workdir", help="toplevel path to where you keep code repos", default="~/code"
    )
    p.add_argument(
        "-wipe", help="wipe before completely recompiling libs", action="store_true"
    )
    p.add_argument(
        "-b", "--buildsys", help="build system (meson or cmake)", default="meson"
    )
    P = p.parse_args()

    dirs = {
        "prefix": Path(P.prefix).expanduser().resolve(),
        "workdir": Path(P.workdir).expanduser().resolve(),
    }

    mumps(P.wipe, dirs, P.buildsys)
