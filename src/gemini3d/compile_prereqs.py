import typing as T
import sys
import os
import pkg_resources
import subprocess
import shutil
import logging
import argparse
from pathlib import Path

from .utils import get_cpu_count
from .web import url_retrieve, extract_tar

BUILDDIR = "build"
NJOBS = get_cpu_count()
# Library parameters

NETCDF_C = "4.7.4"
NETCDF_C_GIT = "https://github.com/Unidata/netcdf-c.git"
NETCDF_FORTRAN = "4.5.2"
NETCDF_FORTRAN_GIT = "https://github.com/Unidata/netcdf-fortran.git"

HDF5_TAG = "1.12/master"
HDF5_GIT = "https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git"

MPIVERSION = "3.1.6"  # OpenMPI 4 needs Scalapack 2.1
MPISHA1 = "bc4cd7fa0a7993d0ae05ead839e6056207e432d4"

LAPACK_GIT = "https://github.com/scivision/lapack.git"
LAPACK_DIR = "lapack"

SCALAPACK_GIT = "https://github.com/scivision/scalapack.git"
SCALAPACK_DIR = "scalapack"

MUMPS_GIT = "https://github.com/scivision/mumps.git"
MUMPS_DIR = "mumps"

# ========= end of user parameters ================

nice = ["nice"] if sys.platform == "linux" else []


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "compiler", help="compiler to build libraries for", choices=["gcc", "intel", "ibmxl"]
    )
    p.add_argument(
        "libs",
        help="libraries to compile",
        choices=["netcdf", "openmpi", "hdf5", "lapack", "scalapack", "mumps"],
        nargs="+",
    )
    p.add_argument("-prefix", help="toplevel path to install libraries under", default="~/lib_gcc")
    p.add_argument("-workdir", help="toplevel path to where you keep code repos", default="~/code")
    p.add_argument("-wipe", help="wipe before completely recompiling libs", action="store_true")
    P = p.parse_args()

    dirs = {
        "prefix": Path(P.prefix).expanduser().resolve(),
        "workdir": Path(P.workdir).expanduser().resolve(),
    }

    setup_libs(P.libs, dirs, P.compiler, P.wipe)


def setup_libs(libs: T.Sequence[str], dirs: T.Dict[str, Path], compiler: str, wipe: bool):

    if compiler == "gcc":
        env = gcc_compilers()
    elif compiler == "intel":
        env = intel_compilers()
    elif compiler == "ibmxl":
        env = ibmxl_compilers()
    else:
        raise ValueError(f"unknown compiler {compiler}")
    # Note: HDF5 needs to be before NetCDF
    if "hdf5" in libs:
        hdf5(dirs, env=env)
    if "netcdf" in libs:
        netcdf_c(dirs, env=env, wipe=wipe)
        netcdf_fortran(dirs, env=env, wipe=wipe)

    # Note: OpenMPI needs to be before scalapack and mumps
    if "openmpi" in libs:
        openmpi(dirs, env=env)
    if "lapack" in libs:
        lapack(wipe, dirs, env=env)
    if "scalapack" in libs:
        scalapack(wipe, dirs, env=env)
    if "mumps" in libs:
        mumps(wipe, dirs, env=env)


def netcdf_c(dirs: T.Dict[str, Path], env: T.Mapping[str, str], wipe: bool = False):
    """ build and install NetCDF-C
    """

    install_dir = dirs["prefix"] / f"netcdf-{NETCDF_C}"
    source_dir = dirs["workdir"] / "netcdf-c"
    build_dir = source_dir / BUILDDIR

    git_update(source_dir, NETCDF_C_GIT, f"v{NETCDF_C}")

    c_args = [
        f"-DCMAKE_INSTALL_PREFIX:PATH={install_dir}",
        "-DCMAKE_BUILD_TYPE:STRING=Release",
        "-DBUILD_SHARED_LIBS:BOOL=ON",
        "-DENABLE_PARALLEL4:BOOL=OFF",
        "-DENABLE_PNETCDF:BOOL=OFF",
        "-DBUILD_UTILITIES:BOOL=OFF",
        "-DENABLE_TESTS:BOOL=off",
        "-DBUILD_TESTING:BOOL=OFF",
        "-DENABLE_HDF4:BOOL=OFF",
        "-DUSE_DAP:BOOL=off",
        "-DENABLE_DAP:BOOL=OFF",
        "-DENABLE_DAP2:BOOL=OFF",
        "-DENABLE_DAP4:BOOL=OFF",
    ]
    cmake_build(c_args, source_dir, build_dir, wipe, env=env, run_test=False)


def netcdf_fortran(dirs: T.Dict[str, Path], env: T.Mapping[str, str], wipe: bool = False):
    """ build and install NetCDF-Fortran
    """

    install_dir = dirs["prefix"] / f"netcdf-{NETCDF_C}"
    source_dir = dirs["workdir"] / "netcdf-fortran"
    build_dir = source_dir / BUILDDIR

    git_update(source_dir, NETCDF_FORTRAN_GIT, f"v{NETCDF_FORTRAN}")

    # NetCDF-Fortran does not yet use NetCDF_ROO
    if sys.platform == "linux":
        netcdf_c = install_dir / "lib/libnetcdf.so"
    elif sys.platform == "win32":
        print(
            "NetCDF4 on MSYS2 may not work, see https://github.com/Unidata/netcdf-c/issues/554",
            file=sys.stderr,
        )
        netcdf_c = install_dir / "bin/libnetcdf.dll"
    elif sys.platform == "darwin":
        netcdf_c = install_dir / "lib/libnetcdf.dylib"
    else:
        raise NotImplementedError(
            f"please open a GitHub Issue for your operating system {sys.platform}"
        )

    patch = [
        f"-DNETCDF_C_LIBRARY:FILEPATH={netcdf_c}",
        f"-DNETCDF_INCLUDE_DIR:PATH={install_dir / 'include'}",
    ]
    f_args = patch + [
        f"-DNetCDF_ROOT:PATH={install_dir}",
        f"-DCMAKE_INSTALL_PREFIX:PATH={install_dir}",
        "-DCMAKE_BUILD_TYPE:STRING=Release",
        "-DBUILD_SHARED_LIBS:BOOL=ON",
        "-DBUILD_UTILITIES:BOOL=OFF",
        "-DENABLE_TESTS:BOOL=off",
        "-DBUILD_EXAMPLES:BOOL=OFF",
    ]
    cmake_build(f_args, source_dir, build_dir, wipe, env=env, run_test=False)


def hdf5(dirs: T.Dict[str, Path], env: T.Mapping[str, str]):
    """ build and install HDF5
    some systems have broken libz and so have trouble extracting tar.bz2 from Python.
    To avoid this, we git clone the release instead.
    """

    if os.name == "nt":
        if "ifort" in env["FC"]:
            msg = """
For Windows with Intel compiler, use HDF5 binaries from HDF Group.
https://www.hdfgroup.org/downloads/hdf5/
look for filename like hdf5-1.12.0-Std-win10_64-vs14-Intel.zip
            """
        elif "gfortran" in env["FC"]:
            msg = """
For MSYS2 on Windows, just use MSYS2 HDF5.
Install from the MSYS2 terminal like:
pacman -S mingw-w64-x86_64-hdf5
reference: https://packages.msys2.org/package/mingw-w64-x86_64-hdf5
            """
        else:
            msg = """
For Windows, use HDF5 binaries from HDF Group.
https://www.hdfgroup.org/downloads/hdf5/
Instead of this, it is generally best to use MSYS2 or Windows Subsystem for Linux
            """
        raise SystemExit(msg)

    hdf5_name = "hdf5"
    install_dir = dirs["prefix"] / hdf5_name
    source_dir = dirs["workdir"] / hdf5_name

    git_update(source_dir, HDF5_GIT, tag=HDF5_TAG)

    cmd = [
        "./configure",
        f"--prefix={install_dir}",
        "--enable-fortran",
        "--enable-build-mode=production",
    ]

    subprocess.check_call(nice + cmd, cwd=source_dir, env=env)

    cmd = ["make", "-C", str(source_dir), "-j", str(NJOBS), "install"]
    subprocess.check_call(nice + cmd)


def openmpi(dirs: T.Dict[str, Path], env: T.Mapping[str, str]):
    """ build and install OpenMPI """
    if os.name == "nt":
        raise NotImplementedError(
            "OpenMPI is not available in native Windows. Use MS-MPI instead."
        )

    mpi_dir = f"openmpi-{MPIVERSION}"
    install_dir = dirs["prefix"] / mpi_dir
    source_dir = dirs["workdir"] / mpi_dir

    tar_name = f"openmpi-{MPIVERSION}.tar.bz2"
    tarfn = dirs["workdir"] / tar_name
    url = f"https://download.open-mpi.org/release/open-mpi/v{MPIVERSION[:3]}/{tar_name}"
    url_retrieve(url, tarfn, ("sha1", MPISHA1))
    extract_tar(tarfn, source_dir)

    cmd = [
        "./configure",
        f"--prefix={install_dir}",
        f"CC={env['CC']}",
        f"CXX={env['CXX']}",
        f"FC={env['FC']}",
    ]

    subprocess.check_call(nice + cmd, cwd=source_dir, env=env)

    cmd = ["make", "-C", str(source_dir), "-j", str(NJOBS), "install"]
    subprocess.check_call(nice + cmd)


def lapack(wipe: bool, dirs: T.Dict[str, Path], env: T.Mapping[str, str]):
    """ build and insall Lapack """
    install_dir = dirs["prefix"] / LAPACK_DIR
    source_dir = dirs["workdir"] / LAPACK_DIR
    build_dir = source_dir / BUILDDIR

    git_update(source_dir, LAPACK_GIT)

    args = [f"-DCMAKE_INSTALL_PREFIX:PATH={install_dir}"]
    cmake_build(args, source_dir, build_dir, wipe, env=env)


def scalapack(wipe: bool, dirs: T.Dict[str, Path], env: T.Mapping[str, str]):
    """ build and install Scalapack """
    source_dir = dirs["workdir"] / SCALAPACK_DIR
    build_dir = source_dir / BUILDDIR

    git_update(source_dir, SCALAPACK_GIT)

    lib_args = [f'-DLAPACK_ROOT={dirs["prefix"] / LAPACK_DIR}']

    args = [f"-DCMAKE_INSTALL_PREFIX:PATH={dirs['prefix'] / SCALAPACK_DIR}"]
    cmake_build(args + lib_args, source_dir, build_dir, wipe, env=env)


def mumps(wipe: bool, dirs: T.Dict[str, Path], env: T.Mapping[str, str]):
    """ build and install Mumps """
    install_dir = dirs["prefix"] / MUMPS_DIR
    source_dir = dirs["workdir"] / MUMPS_DIR
    build_dir = source_dir / BUILDDIR

    scalapack_lib = dirs["prefix"] / SCALAPACK_DIR
    lapack_lib = dirs["prefix"] / LAPACK_DIR

    git_update(source_dir, MUMPS_GIT)

    if env["FC"] == "ifort":
        lib_args = []
    else:
        lib_args = [f"-DSCALAPACK_ROOT:PATH={scalapack_lib}", f"-DLAPACK_ROOT:PATH={lapack_lib}"]

    args = [f"-DCMAKE_INSTALL_PREFIX:PATH={install_dir}"]
    cmake_build(args + lib_args, source_dir, build_dir, wipe, env=env)


def cmake_build(
    args: T.List[str],
    source_dir: Path,
    build_dir: Path,
    wipe: bool,
    env: T.Mapping[str, str],
    run_test: bool = True,
):
    """ build and install with CMake """
    cmake = cmake_minimum_version("3.13")
    cachefile = build_dir / "CMakeCache.txt"
    if wipe and cachefile.is_file():
        cachefile.unlink()

    subprocess.check_call(
        nice + [cmake] + args + ["-B", str(build_dir), "-S", str(source_dir)], env=env
    )

    subprocess.check_call(
        nice + [cmake, "--build", str(build_dir), "--parallel", "--target", "install"]
    )

    if run_test:
        subprocess.check_call(
            nice + ["ctest", "--parallel", "--output-on-failure"], cwd=str(build_dir)
        )


def meson_build(
    args: T.List[str], source_dir: Path, build_dir: Path, wipe: bool, env: T.Mapping[str, str]
) -> int:
    """ build and install with Meson """
    meson = shutil.which("meson")
    if not meson:
        raise FileNotFoundError("Meson not found.")

    if wipe and (build_dir / "build.ninja").is_file():
        args.append("--wipe")

    subprocess.check_call(
        nice + [meson, "setup"] + args + [str(build_dir), str(source_dir)], env=env
    )

    for op in ("test", "install"):
        ret = subprocess.run(nice + [meson, op, "-C", str(build_dir)])

    return ret.returncode


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
    if pkg_resources.parse_version(cmake_ver) < pkg_resources.parse_version(min_version):
        logging.error(f"CMake {cmake_ver} is less than minimum required {min_version}")

    return cmake


def git_update(path: Path, repo: str, tag: str = None):
    """
    Use Git to update a local repo, or clone it if not already existing.

    we use cwd= instead of "git -C" for very old Git versions that might be on your HPC.
    """
    GITEXE = shutil.which("git")

    if not GITEXE:
        logging.error("Git not available.")
        return

    if path.is_dir():
        subprocess.check_call([GITEXE, "-C", str(path), "pull"])
    else:
        # shallow clone
        if tag:
            subprocess.check_call(
                [
                    GITEXE,
                    "clone",
                    repo,
                    "--depth",
                    "1",
                    "--branch",
                    tag,
                    "--single-branch",
                    str(path),
                ]
            )
        else:
            subprocess.check_call([GITEXE, "clone", repo, "--depth", "1", str(path)])


def get_compilers(fc_name: str, cc_name: str, cxx_name: str) -> T.Mapping[str, str]:
    """ get paths to GCC compilers """
    env = os.environ

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


def gcc_compilers() -> T.Mapping[str, str]:
    return get_compilers("gfortran", "gcc", "g++")


def intel_compilers() -> T.Mapping[str, str]:

    fc_name = "ifort"
    cc_name = "icl" if os.name == "nt" else "icc"
    cxx_name = "icl" if os.name == "nt" else "icpc"

    return get_compilers(fc_name, cc_name, cxx_name)


def ibmxl_compilers() -> T.Mapping[str, str]:
    return get_compilers("xlf", "xlc", "xlc++")
