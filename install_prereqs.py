#!/usr/bin/env python
"""
installs libraries
For CentOS HPC systems, consider centos*.sh scripts that use "modules" or "extensions"

For Ubuntu systems, works best for Ubuntu >= 18.04
"""
import subprocess
from pathlib import Path
import sys
from configparser import ConfigParser
from argparse import ArgumentParser
import typing


def os_release() -> typing.List[str]:
    """
    reads /etc/os-release with fallback to legacy methods

    returns
    -------
    'yum' or 'apt'
    """
    fn = Path("/etc/os-release")
    if not fn.is_file():
        if (
            Path("/etc/redhat-release").is_file()
            or Path("/etc/centos-release").is_file()
        ):
            return ["rhel"]
        elif Path("/etc/debian_version").is_file():
            return ["debian"]

    C = ConfigParser(inline_comment_prefixes=("#", ";"))
    ini = "[all]" + fn.read_text()
    C.read_string(ini)
    return C["all"].get("ID_LIKE").strip('"').strip("'").split()


def get_package_manager(like: typing.List[str] = None) -> str:
    if not like:
        like = os_release()
    if isinstance(like, str):
        like = [like]

    if {"centos", "rhel", "fedora"}.intersection(like):
        return "yum"
    elif {"debian", "ubuntu"}.intersection(like):
        return "apt"
    else:
        raise ValueError(
            f"Unknown ID_LIKE={like}, please file bug report or manually specify package manager"
        )


def main(package_manager: str):

    if sys.platform == "linux":
        if not package_manager:
            package_manager = get_package_manager()

        pkgs = {
            "yum": [
                "epel-release",
                "pkg-config",
                "gcc-gfortran",
                "make",
                "ninja-build",
                "MUMPS-openmpi-devel",
                "lapack-devel",
                "scalapack-openmpi-devel",
                "openmpi-devel",
                "octave",
            ],
            "apt": [
                "pkg-config",
                "gfortran",
                "make",
                "ninja-build",
                "libmumps-dev",
                "liblapack-dev",
                "libscalapack-mpi-dev",
                "libopenmpi-dev",
                "octave",
            ],
        }

        if package_manager == "yum":
            if subprocess.run(
                ["sudo", "yum", "--assumeyes", "install"] + pkgs["yum"]
            ).returncode:
                raise SystemExit(
                    "This script is made for personal laptops/desktops.\n"
                    "HPCs using CentOS have system-specific library setup. \n"
                    "If using gfortran, version >= 6 is required."
                    "Try devtoolset-7 if gcc/gfortran is too old on your system."
                )
        elif package_manager == "apt":
            if subprocess.run(["sudo", "apt", "update"]).returncode:
                raise SystemExit("installing prereqs failed.")
            if subprocess.run(
                ["sudo", "apt", "--yes", "install"] + pkgs["apt"]
            ).returncode:
                raise SystemExit("installing prereqs failed.")
    elif sys.platform == "darwin":
        pkgs = {
            "brew": [
                "cmake",
                "gcc",
                "make",
                "ninja",
                "lapack",
                "scalapack",
                "openmpi",
                "octave",
            ]
        }
        if subprocess.run(["brew", "install"] + pkgs["brew"]).returncode:
            raise SystemExit("installing prereqs failed.")
        if subprocess.run(["brew", "tap", "dpo/openblas"]).returncode:
            raise SystemExit("installing prereqs failed.")
        if subprocess.run(["brew", "install", "mumps"]).returncode:
            raise SystemExit("installing prereqs failed.")
    elif sys.platform == "cygwin":
        pkgs = [
            "gcc-fortran",
            "make",
            "ninja",
            "liblapack-devel",
            "libopenmpi-devel",
            "octave",
        ]
        if subprocess.run(["setup-x86_64.exe", "-P"] + pkgs).returncode:
            raise SystemExit("installing prereqs failed.")

        print("use ./compile_prereqs.sh to get Scalapack and Mumps for Cygwin")
    elif sys.platform == "win32":
        raise SystemExit(
            "It is easiest to either use Intel compilers for Windows, or use Windows Subsystem for Linux."
        )
    else:
        raise NotImplementedError(f"unknown platform {sys.platform}")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument(
        "package_manager", help="specify package manager e.g. apt, yum", nargs="?"
    )
    P = p.parse_args()

    main(P.package_manager)
