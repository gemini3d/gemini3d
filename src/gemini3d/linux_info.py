#!/usr/bin/env python
"""
get Linux system info
"""
from configparser import ConfigParser
import typing as T
from pathlib import Path
import sys


def os_release() -> T.List[str]:
    """
    reads /etc/os-release with fallback to legacy methods

    Returns
    -------

    linux_names: list of str
        name(s) of operating system detect
    """

    if sys.platform != "linux":
        return []

    fn = Path("/etc/os-release")
    if not fn.is_file():
        if Path("/etc/redhat-release").is_file() or Path("/etc/centos-release").is_file():
            return ["rhel"]
        elif Path("/etc/debian_version").is_file():
            return ["debian"]

    return parse_os_release("[all]" + fn.read_text())


def parse_os_release(txt: str) -> T.List[str]:
    """ parse /etc/os-release text """

    C = ConfigParser(inline_comment_prefixes=("#", ";"))
    C.read_string(txt)
    return C["all"].get("ID_LIKE").strip('"').strip("'").split()


def get_package_manager(like: T.List[str] = None) -> str:
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
