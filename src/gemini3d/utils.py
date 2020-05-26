import subprocess
import os
import shutil
from pathlib import Path
from datetime import datetime, timedelta
import typing as T

try:
    import psutil
except ImportError:
    psutil = None
    # pip install psutil will improve CPU utilization.


git = shutil.which("git")

Pathlike = T.Union[str, Path]

__all__ = ["gitrev", "get_cpu_count", "ymdhourdec2datetime", "datetime2ymd_hourdec"]


def gitrev() -> str:
    if not git:
        return ""

    return subprocess.check_output(
        [git, "rev-parse", "--short", "HEAD"], universal_newlines=True
    ).strip()


def get_cpu_count() -> int:
    """ get a physical CPU count

    Returns
    -------
    count: int
        detect number of physical CPU
    """

    max_cpu = None
    # without psutil, hyperthreaded CPU may overestimate physical count by factor of 2 (or more)
    if psutil is not None:
        max_cpu = psutil.cpu_count(logical=False)
        extradiv = 1
        if max_cpu is None:
            max_cpu = psutil.cpu_count()
            extradiv = 2
    if max_cpu is None:
        max_cpu = os.cpu_count()
        extradiv = 2

    return max_cpu // extradiv


def ymdhourdec2datetime(year: int, month: int, day: int, hourdec: float) -> datetime:
    """
    convert year,month,day + decimal hour HH.hhh to time
    """

    return datetime(year, month, day, int(hourdec), int((hourdec * 60) % 60)) + timedelta(
        seconds=(hourdec * 3600) % 60
    )


def datetime2ymd_hourdec(dt: datetime) -> str:
    """
    convert datetime to ymd_hourdec string for filename stem
    """

    return (
        dt.strftime("%Y%m%d")
        + f"_{dt.hour*3600 + dt.minute*60 + dt.second + dt.microsecond/1e6:12.6f}"
    )
