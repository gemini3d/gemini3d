#!/usr/bin/env python3
"""
Installs prereqs for Gemini program for Gfortran
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
from meson_file_download import url_retrieve
from meson_file_extract import extract_tar
from meson_cpu_count import get_cpu_count

# ========= user parameters ======================
BUILDDIR = "build"
NJOBS = get_cpu_count()
# Library parameters
HDF5VERSION = "1.10.6"
HDF5URL = "https://zenodo.org/record/3659270/files/hdf5-1.10.6.tar.bz2?download=1"
HDF5MD5 = "03095102a6118c32a75a9b9b40be66f2"
HDF5DIR = f"hdf5-{HDF5VERSION}"

MPIVERSION = "3.1.5"  # OpenMPI 4 needs Scalapack 2.1
MPISHA1 = "56a74b116c81d4f3704c051a67e4422094ff913d"
MPIDIR = f"openmpi-{MPIVERSION}"

LAPACKGIT = "https://github.com/scivision/lapack"
LAPACKDIR = "lapack"

SCALAPACKGIT = "https://github.com/scivision/scalapack"
SCALAPACKDIR = "scalapack"

MUMPSGIT = "https://github.com/scivision/mumps"
MUMPSDIR = "mumps"

# ========= end of user parameters ================

nice = ["nice"] if sys.platform == "linux" else []


def hdf5(dirs: typing.Dict[str, Path]):
    if os.name == "nt":
        raise SystemExit("Please use binaries from HDF Group for Windows appropriate for your compiler.")

    install_dir = dirs["prefix"] / HDF5DIR
    source_dir = dirs["workdir"] / HDF5DIR

    tarfn = f"hdf5-{HDF5VERSION}.tar.bz2"
    url_retrieve(HDF5URL, tarfn, ("md5", HDF5MD5))
    extract_tar(tarfn, source_dir)

    env = get_compilers()

    cmd = nice + ["./configure", f"--prefix={install_dir}", "--enable-fortran", "--enable-build-mode=production"]

    subprocess.check_call(cmd, cwd=source_dir, env=env)

    cmd = nice + ["make", "-C", str(source_dir), f"-j{NJOBS}", "install"]
    subprocess.check_call(cmd)


def openmpi(dirs: typing.Dict[str, Path]):
    if os.name == "nt":
        raise SystemExit("OpenMPI is not available in native Windows.")

    install_dir = dirs["prefix"] / MPIDIR
    source_dir = dirs["workdir"] / MPIDIR

    tarfn = f"openmpi-{MPIVERSION}.tar.bz2"
    url = f"https://download.open-mpi.org/release/open-mpi/v{MPIVERSION[:3]}/{tarfn}"
    url_retrieve(url, tarfn, ("sha1", MPISHA1))
    extract_tar(tarfn, source_dir)

    env = get_compilers()

    cmd = nice + ["./configure", f"--prefix={install_dir}", f"CC={env['CC']}", f"CXX={env['CXX']}", f"FC={env['FC']}"]

    subprocess.check_call(cmd, cwd=source_dir, env=env)

    cmd = nice + ["make", "-C", str(source_dir), f"-j{NJOBS}", "install"]
    subprocess.check_call(cmd)


def lapack(wipe: bool, dirs: typing.Dict[str, Path], buildsys: str):
    install_dir = dirs["prefix"] / LAPACKDIR
    source_dir = dirs["workdir"] / LAPACKDIR
    build_dir = source_dir / BUILDDIR

    update(source_dir, LAPACKGIT)

    if buildsys == "cmake":
        args = [f"-DCMAKE_INSTALL_PREFIX={install_dir}"]
        cmake_build(args, source_dir, build_dir, wipe, env=get_compilers())
    elif buildsys == "meson":
        args = [f"--prefix={dirs['prefix']}"]
        meson_build(args, source_dir, build_dir, wipe, env=get_compilers())
    else:
        raise ValueError(f"unknown build system {buildsys}")


def scalapack(wipe: bool, dirs: typing.Dict[str, Path], buildsys: str):
    source_dir = dirs["workdir"] / SCALAPACKDIR
    build_dir = source_dir / BUILDDIR

    update(source_dir, SCALAPACKGIT)

    lib_args = [f'-DLAPACK_ROOT={dirs["prefix"] / LAPACKDIR}']

    if buildsys == "cmake":
        args = [f"-DCMAKE_INSTALL_PREFIX={dirs['prefix'] / SCALAPACKDIR}"]
        cmake_build(args + lib_args, source_dir, build_dir, wipe, env=get_compilers())
    elif buildsys == "meson":
        args = [f"--prefix={dirs['prefix']}"]
        meson_build(args + lib_args, source_dir, build_dir, wipe, env=get_compilers())

    else:
        raise ValueError(f"unknown build system {buildsys}")


def mumps(wipe: bool, dirs: typing.Dict[str, Path], buildsys: str):
    install_dir = dirs["prefix"] / MUMPSDIR
    source_dir = dirs["workdir"] / MUMPSDIR
    build_dir = source_dir / BUILDDIR

    scalapack_lib = dirs["prefix"] / SCALAPACKDIR
    lapack_lib = dirs["prefix"] / LAPACKDIR

    update(source_dir, MUMPSGIT)

    lib_args = [f"-DSCALAPACK_ROOT={scalapack_lib}", f"-DLAPACK_ROOT={lapack_lib}"]

    if buildsys == "cmake":
        args = [f"-DCMAKE_INSTALL_PREFIX={install_dir}"]
        cmake_build(args + lib_args, source_dir, build_dir, wipe, env=get_compilers())
    elif buildsys == "meson":
        args = [f"--prefix={dirs['prefix']}"]
        meson_build(args + lib_args, source_dir, build_dir, wipe, env=get_compilers())
    else:
        raise ValueError(f"unknown build system {buildsys}")


def cmake_build(args: typing.List[str], source_dir: Path, build_dir: Path, wipe: bool, env: typing.Mapping[str, str]):
    cmake = cmake_minimum_version("3.13")
    cachefile = build_dir / "CMakeCache.txt"
    if wipe and cachefile.is_file():
        cachefile.unlink()

    subprocess.check_call(nice + [cmake] + args + ["-B", str(build_dir), "-S", str(source_dir)], env=env)

    subprocess.check_call(nice + [cmake, "--build", str(build_dir), "--parallel", "--target", "install"])

    subprocess.check_call(nice + ["ctest", "--parallel", "--output-on-failure"], cwd=str(build_dir))


def meson_build(args: typing.List[str], source_dir: Path, build_dir: Path, wipe: bool, env: typing.Mapping[str, str]):
    meson = shutil.which("meson")
    if not meson:
        raise FileNotFoundError("Meson not found.")

    if wipe and (build_dir / "build.ninja").is_file():
        args.append("--wipe")

    subprocess.check_call(nice + [meson, "setup"] + args + [str(build_dir), str(source_dir)], env=env)

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
    cmake_ver = subprocess.check_output([cmake, "--version"], universal_newlines=True).split("\n")[0].split(" ")[2]
    if pkg_resources.parse_version(cmake_ver) < pkg_resources.parse_version(min_version):
        logging.error(f"CMake {cmake_ver} is less than minimum required {min_version}")

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
        subprocess.run([GITEXE, "pull"], cwd=str(path))
    else:
        subprocess.run([GITEXE, "clone", repo, str(path)])


@lru_cache()
def get_compilers() -> typing.Mapping[str, str]:

    env = os.environ

    fc_name = "gfortran"
    cc_name = "gcc"
    cxx_name = "g++"

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
    p.add_argument("libs", help="libraries to compile", choices=["openmpi", "hdf5", "lapack", "scalapack", "mumps"], nargs="+")
    p.add_argument("-prefix", help="toplevel path to install libraries under", default="~/lib_gcc")
    p.add_argument("-workdir", help="toplevel path to where you keep code repos", default="~/code")
    p.add_argument("-wipe", help="wipe before completely recompiling libs", action="store_true")
    p.add_argument("-b", "--buildsys", help="build system (meson or cmake)", default="cmake")
    P = p.parse_args()

    dirs = {"prefix": Path(P.prefix).expanduser().resolve(), "workdir": Path(P.workdir).expanduser().resolve()}

    if "hdf5" in P.libs:
        hdf5(dirs)
    if "openmpi" in P.libs:
        openmpi(dirs)
    if "lapack" in P.libs:
        lapack(P.wipe, dirs, P.buildsys)
    if "scalapack" in P.libs:
        scalapack(P.wipe, dirs, P.buildsys)
    if "mumps" in P.libs:
        mumps(P.wipe, dirs, P.buildsys)
