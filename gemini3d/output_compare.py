"""
compare simulation outputs to verify model performance
"""
from pathlib import Path
import numpy as np
import logging
import sys
import typing

from .readdata import read_config, loadframe, datetime_range

try:
    from .plotdiff import plotdiff
except ModuleNotFoundError:
    plotdiff = None


def compare_all(outdir: Path, refdir: Path, tol: typing.Dict[str, float], doplot: bool = False) -> int:
    """
    compare two directories across time steps
    """
    outdir = Path(outdir).expanduser()
    refdir = Path(refdir).expanduser()

    if not outdir.is_dir():
        print(f"output directory(s) not found: {outdir}", file=sys.stderr)
        raise SystemExit(77)
    if not refdir.is_dir():
        print(f"reference directory(s) not found: {refdir}", file=sys.stderr)
        raise SystemExit(77)
    if outdir.samefile(refdir):
        raise OSError(f"reference and output are the same directory: {outdir}")

    ref: typing.Dict[str, typing.Any] = {}
    # %% READ IN THE SIMULATION INFORMATION
    outinput_dir = outdir / "inputs"
    params = read_config(outinput_dir)
    # params["lxs"] = get_simsize(outinput_dir)
    # %% TIMES OF INTEREST
    t0 = params["t0"]
    times = datetime_range(t0, t0 + params["tdur"], params["dtout"])
    if len(times) <= 1:
        raise ValueError(f"{outdir} simulation did not run long enough, must run for more than one time step")

    errs = 0
    for i, t in enumerate(times):
        st = f"UTsec {t}"
        A = loadframe(outdir, t)
        B = loadframe(refdir, t)

        names = ["ne", "v1", "v2", "v3", "Ti", "Te", "J1", "J2", "J3"]
        itols = ["N", "V", "V", "V", "T", "T", "J", "J", "J"]

        for k, j in zip(names, itols):
            a = A[k][1]
            b = B[k][1]
            assert a.shape == b.shape, f"{k} time {i} {t}: ref shape {b.shape} does not match data shape {a.shape}"
            if not np.allclose(a, b, tol[f"rtol{j}"], tol[f"atol{j}"], True):
                errs += 1
                logging.error(f"{k} {st}   {abs(a - b).max().item():.3e}")
                if plotdiff is not None:
                    plotdiff(a, b, k, t, outdir, refdir)
        # %% assert time steps have unique output (earth always rotating...)
        if i > 1:
            names = ["ne", "v1", "v2", "v3"]
            itols = ["N", "V", "V", "V"]
            for k, j in zip(names, itols):
                if np.allclose(ref[k][1], a, 0.0001 * tol[f"rtol{j}"], 0.0001 * tol[f"atol{j}"]):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        if i == 3:
            for k in ("Ti", "Te"):
                if np.allclose(ref[k][1], A[k][1], 0.01 * tol["rtolT"], 0.1 * tol["atolT"]):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        if i == 2:
            for k in ("J1", "J2", "J3"):
                if np.allclose(ref[k][1], a, 0.01 * tol["rtolJ"], 0.1 * tol["atolJ"]):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        ref.update(A)

    return errs
