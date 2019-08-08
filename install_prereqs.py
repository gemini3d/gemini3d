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
    'rhel' or 'debian'
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
                "ninja-build",
                "MUMPS-openmpi-devel",
                "lapack-devel",
                "blacs-openmpi-devel",
                "scalapack-openmpi-devel",
                "openmpi-devel",
            ],
            "apt": [
                "pkg-config",
                "gfortran",
                "ninja-build",
                "libmumps-dev",
                "liblapack-dev",
                "libblacs-mpi-dev",
                "libscalapack-mpi-dev",
                "libopenmpi-dev",
                "openmpi-bin",
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
        else:
            raise ValueError(
                f"I don't know package manager {package_manager}, try installing the prereqs manually"
            )
    elif sys.platform == "darwin":
        pkgs = {"brew": ["gcc", "ninja", "lapack", "scalapack", "openmpi"]}
        if subprocess.run(["brew", "install"] + pkgs["brew"]).returncode:
            raise SystemExit("installing prereqs failed.")
        # if subprocess.run(["brew", "tap", "dpo/openblas"]).returncode:
        #    raise SystemExit("installing prereqs failed.")
        # if subprocess.run(["brew", "install", "mumps"]).returncode:
        #    raise SystemExit("installing prereqs failed.")
    elif sys.platform == "cygwin":
        pkgs = ["gcc-fortran", "ninja", "liblapack-devel", "libopenmpi-devel"]
        if subprocess.run(["setup-x86_64.exe", "-P"] + pkgs).returncode:
            raise SystemExit("installing prereqs failed.")

        print("meson will automatically get Scalapack and Mumps during Gemini build")
    elif sys.platform == "win32":
        raise SystemExit(
            "It is easiest to use Intel compilers for Windows, or Windows Subsystem for Linux."
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

    print('If you need "meson" do "python3 -m pip install --user meson"')
