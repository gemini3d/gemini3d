#!/usr/bin/env python3
"""
installs Gemini prerequisite libraries for CentOS, Debian, Ubuntu, Homebrew and Cygwin, assuming GCC/Gfortran

Michael Hirsch, Ph.D.
"""
import subprocess
import sys
from argparse import ArgumentParser
import typing as T

from gemini.linux_info import get_package_manager


def main(package_manager: str):

    cmd: T.List[str] = []  # so that it's sure to be defined

    if sys.platform == "linux":
        if not package_manager:
            package_manager = get_package_manager()

        pkgs = {
            "yum": [
                "epel-release",
                "pkg-config",
                "gcc-gfortran",
                "MUMPS-openmpi-devel",
                "lapack-devel",
                "blacs-openmpi-devel",
                "scalapack-openmpi-devel",
                "openmpi-devel",
            ],
            "apt": [
                "pkg-config",
                "gfortran",
                "libmumps-dev",
                "liblapack-dev",
                "libblacs-mpi-dev",
                "libscalapack-mpi-dev",
                "libopenmpi-dev",
                "openmpi-bin",
            ],
        }

        if package_manager == "yum":
            subprocess.run(["sudo", "yum", "--assumeyes", "updateinfo"])
            cmd = ["sudo", "yum", "--assumeyes", "install"] + pkgs["yum"]
        elif package_manager == "apt":
            subprocess.run(["sudo", "apt", "update", "--yes"])
            cmd = ["sudo", "apt", "--yes", "install"] + pkgs["apt"]
        else:
            raise ValueError(f"Unknown package manager {package_manager}, try installing the prereqs manually")
    elif sys.platform == "darwin":
        pkgs = {"brew": ["gcc", "make", "cmake", "lapack", "openmpi"]}
        cmd = ["brew", "install"] + pkgs["brew"]
        # just autobuild Mumps instead, it's much faster
        # subprocess.run(["brew", "tap", "dpo/openblas"])
        # subprocess.run(["brew", "install", "mumps"])
    elif sys.platform == "cygwin":
        pkgs = ["gcc-fortran", "liblapack-devel", "libopenmpi-devel"]
        cmd = ["setup-x86_64.exe", "-P"] + pkgs
    elif sys.platform == "win32":
        raise SystemExit("Windows Subsystem for Linux is recommended.")
    else:
        raise ValueError(f"unknown platform {sys.platform}")

    print(cmd)
    ret = subprocess.run(cmd)

    raise SystemExit(ret.returncode)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("package_manager", help="specify package manager e.g. apt, yum", nargs="?")
    P = p.parse_args()

    main(P.package_manager)
