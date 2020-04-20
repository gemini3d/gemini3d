import os
import shutil
import subprocess
import functools
import typing as T
from pathlib import Path

Pathlike = T.Union[str, Path]


@functools.lru_cache()
def wsl_available() -> bool:
    """
    heuristic to detect if Windows Subsystem for Linux is available.

    Uses presence of /etc/os-release in the WSL image to say Linux is there.
    This is a de facto file standard across Linux distros.
    """

    has_wsl = False
    if os.name == "nt" and shutil.which("wsl"):
        has_wsl = wsl_file_exist("/etc/os-release")

    return has_wsl


def wsl_file_exist(file: Pathlike) -> bool:
    """
    path is specified as if in WSL
    NOT //wsl$/Ubuntu/etc/os-release
    but /etc/os-release
    """
    if os.name != "nt":
        return False

    try:
        return subprocess.run(["wsl", "test", "-f", str(file)], timeout=10).returncode == 0
    except subprocess.TimeoutExpired:
        return False
