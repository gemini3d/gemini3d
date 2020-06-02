from matplotlib.figure import Figure
import numpy as np
from pathlib import Path
from datetime import datetime

from .utils import gitrev


def plotdiff(A: np.ndarray, B: np.ndarray, name: str, time: datetime, outdir: Path, refdir: Path):

    A = A.squeeze()
    B = B.squeeze()

    if A.ndim == 3:
        if A.shape[0] == 7:
            # presuming 7 => species and 7th (last) index is electrons
            A = A[-1, :, :]
            B = B[-1, :, :]

    fg = Figure(tight_layout=True, figsize=(12, 5))
    axs = fg.subplots(1, 3)

    if A.ndim == 2:
        diff2d(A, B, name, fg, axs)
    elif A.ndim == 1:
        diff1d(A, B, name, fg, axs)
    else:
        fg.close()
        print(f"skipping diff plot: {name}")
        return None

    axs[0].set_title(str(outdir))
    axs[1].set_title(str(refdir))
    axs[2].set_title(f"diff: {name}")

    tstr = time.isoformat()
    ttxt = f"{name}  {tstr}  Git: {gitrev()}"

    fg.suptitle(ttxt)

    fn = outdir / f"{name}-diff-{tstr.replace(':','')}.png"
    print("writing", fn)
    fg.savefig(fn)


def diff1d(A: np.ndarray, B: np.ndarray, name: str, fg, axs):

    axs[0].plot(A)

    axs[1].plot(B)

    axs[2].plot(A - B)


def diff2d(A: np.ndarray, B: np.ndarray, name: str, fg, axs):

    cmap = "bwr" if name.startswith(("J", "v")) else None

    hi = axs[0].pcolormesh(A, cmap=cmap)
    fg.colorbar(hi, ax=axs[0])

    hi = axs[1].pcolormesh(B, cmap=cmap)
    fg.colorbar(hi, ax=axs[1])

    hi = axs[2].pcolormesh(A - B, cmap="bwr")
    fg.colorbar(hi, ax=axs[2])
