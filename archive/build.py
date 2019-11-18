#!/usr/bin/env python3
"""
build.py makes building projects with CMake or Meson + Ninja even simpler.
It facilitates easy testing across operating systems and compiler vendors.
Michael Hirsch, Ph.D.

## Per-compiler tips

### PGI

PATH must include the PGI compilers bin/ directory before running build.py.

### Intel

The Intel compiler environment must be configured before running build.py:

* Windows: compilervars.bat intel64
* Linux / Mac: source compilervars.sh intel64
"""
from pathlib import Path
import os
import sys
import shutil
import subprocess
from typing import Dict, List, Tuple
from argparse import ArgumentParser

if sys.version_info < (3, 6):
    raise RuntimeError("build.py requires Python >= 3.6")

MESON = shutil.which("meson")
NINJA = shutil.which("ninja")
CMAKE = shutil.which("cmake")
CTEST = shutil.which("ctest")

MSVC = "Visual Studio 15 2017"

# Must have .resolve() to work in general regardless of invocation directory
SRC = Path(__file__).parent.resolve()
BUILD = SRC / "build"

LIBPREFIX = "~/lib_"
MUMPSDIR = "mumps"


def do_build(buildsys: str, compilers: Dict[str, str], args: List[str], **kwargs):
    """
    attempts build with Meson or CMake
    """
    if buildsys == "meson" and MESON and NINJA:
        meson_setup(compilers, args, **kwargs)
    elif buildsys == "cmake" and CMAKE:
        cmake_setup(compilers, args, **kwargs)
    else:
        raise FileNotFoundError("Could not find CMake or Meson + Ninja")


def _needs_wipe(fn: Path, wipe: bool) -> bool:
    """
    This detection of regeneration needed is not perfect.
    """
    if not fn.is_file():
        return False

    if wipe:
        return True

    with fn.open() as f:
        for line in f:
            if line.startswith("CMAKE_C_COMPILER:FILEPATH"):
                cc = line.split("/")[-1].strip()  # must have strip() for junk in cache
                if cc != compilers["CC"]:
                    print("regenerating due to C compiler change:", cc, "=>", compilers["CC"])
                    wipe = True
                    break
            elif line.startswith("CMAKE_GENERATOR:INTERNAL"):
                gen = line.split("=")[-1]
                if gen.startswith("Unix") and os.name == "nt":
                    print("regenerating due to OS change: Unix => Windows")
                    wipe = True
                    break
                elif gen.startswith(("MinGW", "Visual")) and os.name != "nt":
                    print("regenerating due to OS change: Windows => Unix")
                    wipe = True
                    break
                elif gen.startswith("Visual") and compilers["CC"] != "cl":
                    print("regenerating due to C compiler change: MSVC =>", compilers["CC"])
                    wipe = True
                    break

    return wipe


def cmake_setup(compilers: Dict[str, str], args: List[str], **kwargs):
    """
    attempt to build using CMake >= 3
    """
    if compilers["CC"] == "cl":
        wopts = ["-G", MSVC, "-A", "x64"]
    elif os.name == "nt":
        wopts = ["-G", "MinGW Makefiles", "-DCMAKE_SH=CMAKE_SH-NOTFOUND"]
    else:
        wopts = []

    wopts += args

    if kwargs.get("install"):  # path specified
        wopts.append("-DCMAKE_INSTALL_PREFIX:PATH=" + str(Path(kwargs["install"]).expanduser()))

    cachefile = BUILD / "CMakeCache.txt"

    if _needs_wipe(cachefile, kwargs.get("wipe")):
        cachefile.unlink()
        shutil.rmtree(BUILD / "CMakeFiles", ignore_errors=True)

    # we didn't use -S -B to be compatible with CMake < 3.12
    ret = subprocess.run([CMAKE] + wopts + [str(SRC)], cwd=BUILD, env=os.environ.update(compilers))
    if ret.returncode:
        raise SystemExit(ret.returncode)

    ret = subprocess.run([CMAKE, "--build", str(BUILD), "--parallel"])

    test_result(ret)

    # %% test
    _cmake_test(kwargs.get("dotest"))
    # %% install
    if kwargs.get("install"):
        subprocess.run([CMAKE, "--build", str(BUILD), "--parallel", "--target", "install"])
        if ret.returncode:
            raise SystemExit(ret.returncode)


def _cmake_test(dotest: bool):
    if not dotest:
        return

    if not CTEST:
        raise FileNotFoundError("CTest not available")

    if compilers["CC"] == "cl":
        ret = subprocess.run([CMAKE, "--build", str(BUILD), "--target", "RUN_TESTS"])
        if ret.returncode:
            raise SystemExit(ret.returncode)
    else:
        ret = subprocess.run([CTEST, "--parallel", "--output-on-failure"], cwd=BUILD)
        if ret.returncode:
            raise SystemExit(ret.returncode)


def meson_setup(compilers: Dict[str, str], args: List[str], **kwargs):
    """
    attempt to build with Meson + Ninja
    """
    build_ninja = BUILD / "build.ninja"

    meson_setup = [MESON] + ["setup"] + args

    if kwargs.get("install"):
        meson_setup.append("--prefix " + str(Path(kwargs["install"]).expanduser()))

    if kwargs.get("wipe") and build_ninja.is_file():
        meson_setup.append("--wipe")

    meson_setup += [str(BUILD), str(SRC)]

    if kwargs.get("wipe") or not build_ninja.is_file():
        ret = subprocess.run(meson_setup, env=os.environ.update(compilers))
        if ret.returncode:
            raise SystemExit(ret.returncode)

    ret = subprocess.run([NINJA, "-C", str(BUILD)])

    test_result(ret)

    if kwargs.get("dotest"):
        if not ret.returncode:
            ret = subprocess.run([MESON, "test", "-C", str(BUILD)])  # type: ignore     # MyPy bug
            if ret.returncode:
                raise SystemExit(ret.returncode)

    if kwargs.get("install"):
        if not ret.returncode:
            ret = subprocess.run([MESON, "install", "-C", str(BUILD)])  # type: ignore     # MyPy bug
            if ret.returncode:
                raise SystemExit(ret.returncode)


def test_result(ret: subprocess.CompletedProcess):
    if not ret.returncode:
        print("\nBuild Complete!")
    else:
        raise SystemExit(ret.returncode)


# %% compilers
def clang_params() -> Tuple[Dict[str, str], List[str]]:
    """
    LLVM compilers e.g. Clang, Flang
    """
    compilers = {"CC": "clang", "CXX": "clang++", "FC": "flang"}

    args: List[str] = []

    return compilers, args


def gnu_params() -> Tuple[Dict[str, str], List[str]]:
    """
    GNU compilers e.g. GCC, Gfortran
    """
    compilers = {"FC": "gfortran", "CC": "gcc", "CXX": "g++"}

    # determine library dir
    libprefix = Path(LIBPREFIX + "gcc").expanduser()

    args: List[str] = [f"-DMUMPS_ROOT={libprefix / MUMPSDIR}"]

    return compilers, args


def intel_params() -> Tuple[Dict[str, str], List[str]]:
    """
    Intel compilers
    """
    if not os.environ.get("MKLROOT"):
        raise EnvironmentError("must have set MKLROOT by running compilervars.bat or source compilervars.sh before this script.")

    # %% compiler variables
    compilers = {"FC": "ifort"}

    if os.name == "nt":
        compilers["CC"] = compilers["CXX"] = "icl.exe"
    else:
        compilers["CC"] = "icc"
        compilers["CXX"] = "icpc"

    # determine library dir
    libprefix = Path(LIBPREFIX + "intel").expanduser()

    args: List[str] = [f"-DMUMPS_ROOT={libprefix / MUMPSDIR}"]

    return compilers, args


def msvc_params() -> Tuple[Dict[str, str], List[str]]:
    """
    Micro$oft Visual Studio

    Note in general MSVC doesn't have good modern C++ features,
    so don't be surprised if a C++11 or newer program doesn't compile.
    """
    if not shutil.which("cl"):
        raise EnvironmentError("Must have PATH set to include MSVC cl.exe compiler bin directory")

    compilers = {"CC": "cl", "CXX": "cl"}

    args: List[str] = []

    return compilers, args


def pgi_params() -> Tuple[Dict[str, str], List[str]]:
    """
    Nvidia PGI compilers

    pgc++ is not available on Windows at this time
    """
    if not shutil.which("pgcc") or not shutil.which("pgfortran"):
        raise EnvironmentError("Must have PATH set to include PGI compiler bin directory")

    # %% compiler variables
    compilers = {"FC": "pgfortran", "CC": "pgcc"}
    if os.name != "nt":
        compilers["CXX"] = "pgc++"

    args: List[str] = []

    return compilers, args


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("vendor", help="compiler vendor [clang, gnu, intel, msvc, pgi]", nargs="?", default="gnu")
    p.add_argument("-wipe", help="wipe and rebuild from scratch", action="store_true")
    p.add_argument("-buildsys", help="default build system", default="cmake")
    p.add_argument("-args", help="preprocessor arguments", nargs="+", default=[])
    p.add_argument("-debug", help="debug (-O0) instead of release (-O3) build", action="store_true")
    p.add_argument("-test", help="run self-test / example", action="store_true")
    p.add_argument("-install", help="specify full install directory e.g. ~/lib_gcc/mumps")
    a = p.parse_args()

    if a.vendor == "clang":
        compilers, args = clang_params()
    elif a.vendor in ("gnu", "gcc"):
        compilers, args = gnu_params()
    elif a.vendor == "intel":
        compilers, args = intel_params()
    elif a.vendor == "msvc":
        compilers, args = msvc_params()
    elif a.vendor == "pgi":
        compilers, args = pgi_params()
    else:
        raise ValueError(a.vendor)

    args += a.args
    if a.debug:
        args.append("-DCMAKE_BUILD_TYPE=Debug")

    do_build(a.buildsys, compilers, args, wipe=a.wipe, dotest=a.test, install=a.install)
