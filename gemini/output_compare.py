"""
compare simulation outputs to verify model performance
"""
from pathlib import Path
import numpy as np
import logging
import sys
import typing

from .readdata import readconfig, loadframe, datetime_range


def compare_all(dir1: Path, dir2: Path, tol: typing.Dict[str, float]) -> int:
    """
    compare two directories across time steps
    """
    dir1 = Path(dir1).expanduser()
    dir2 = Path(dir2).expanduser()

    if not dir1.is_dir() or not dir2.is_dir():
        print("comparison directory(s) not found: {dir1}  {dir2}", file=sys.stderr)
        raise SystemExit(77)

    if dir1.samefile(dir2):
        raise OSError(f"compare_all inputs are the same directory: {dir1}")

    ref: typing.Dict[str, typing.Any] = {}

    # %% READ IN THE SIMULATION INFORMATION
    params = readconfig(dir1 / "inputs/config.ini")
    # %% TIMES OF INTEREST
    t0 = params["t0"]
    times = datetime_range(t0, t0 + params["tdur"], params["dtout"])
    if len(times) <= 1:
        raise ValueError(
            "simulation did not run long enough, must run for more than one time step"
        )

    errs = 0

    for i, t in enumerate(times):
        st = f"UTsec {t}"
        A = loadframe(dir1, t)
        B = loadframe(dir2, t)

        if not np.allclose(A["ne"][1], B["ne"][1], tol["rtolN"], tol["atolN"], True):
            errs += 1
            logging.error(f"Ne {st}   {abs(A['ne'][1] - B['ne'][1]).max().item():.3e}")

        for k in ("v1", "v2", "v3"):
            if not np.allclose(A[k][1], B[k][1], tol["rtolV"], tol["atolV"], True):
                errs += 1
                logging.error(f"{k} {st}   {abs(A[k][1] - B[k][1]).max().item():.3e}")

        for k in ("Ti", "Te"):
            if not np.allclose(A[k][1], B[k][1], tol["rtolT"], tol["atolT"], True):
                errs += 1
                logging.error(f"{k} {st}   {abs(A[k][1] - B[k][1]).max().item():.3e}")

        for k in ("J1", "J2", "J3"):
            if not np.allclose(A[k][1], B[k][1], tol["rtolJ"], tol["atolJ"], True):
                errs += 1
                logging.error(f"{k} {st}  {abs(A[k][1] - B[k][1]).max().item():.3e}")

        # %% assert time steps have unique output (earth always rotating...)
        if i > 1:
            if np.allclose(ref["ne"][1], A["ne"][1], tol["rtol"], tol["atol"]):
                errs += 1
                logging.error(f"Ne {st} too similar to prior step")

            for k in ("v1", "v2", "v3"):
                if np.allclose(ref[k][1], A[k][1], tol["rtol"], tol["atol"]):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        if i == 3:
            for k in ("Ti", "Te"):
                if np.allclose(ref[k][1], A[k][1], tol["rtol"], tol["atol"]):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        if i == 2:
            for k in ("J1", "J2", "J3"):
                if np.allclose(ref[k][1], A[k][1], tol["rtol"], tol["atol"]):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        ref.update(A)

    return errs
