from matplotlib.figure import Figure
import numpy as np
from pathlib import Path
from datetime import datetime
import shutil
import subprocess

git = shutil.which("git")


def plotdiff(A: np.ndarray, B: np.ndarray, name: str, time: datetime, dir1: Path, dir2: Path):
    A = A.squeeze()
    B = B.squeeze()
    if A.ndim != 2 or B.ndim != 2:
        return None

    fg = Figure(tight_layout=True, figsize=(12, 5))
    axs = fg.subplots(1, 2)

    hi = axs[0].pcolormesh(A)
    fg.colorbar(hi, ax=axs[0])
    axs[0].set_title(str(dir1))

    hi = axs[1].pcolormesh(B)
    fg.colorbar(hi, ax=axs[1])
    axs[1].set_title(str(dir2))

    ttxt = f"{name}  {time.isoformat()}"
    if git:
        ttxt += (
            "  Git: "
            + subprocess.check_output(
                [git, "rev-parse", "--short", "HEAD"], universal_newlines=True
            ).strip()
        )

    fg.suptitle(ttxt)

    fn = dir2 / f"{name}-diff-{time.isoformat().replace(':','')}.png"
    print("writing", fn)
    fg.savefig(fn)
