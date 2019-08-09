"""
compare simulation outputs to verify model performance
"""
from pathlib import Path
import numpy as np
import logging
import typing

from .readdata import readconfig, loadframe, datetime_range


def compare_all(dir1: Path, dir2: Path) -> int:
    """
    the absolute and relative tolerance account for slight IEEE-754 based differences,
    including non-associativity of floating-point arithmetic.
    these parameters are a bit arbitrary.

    % per MZ Oct 17, 2018:
    % Ti,Te=1 K
    % ne=1e6 m-3
    % vi,v2,v3=1 m/s
    % J1,J2,J3 = 1e-9

    % MZ wants to change what we consider signficant...
    % Ti,Te=5 K
    % ne=1e7 m-3
    % vi,v2,v3=2 m/s
    % J1,J2,J3 = 1e-9
    """

    dir1 = Path(dir1).expanduser().resolve(strict=True)
    dir2 = Path(dir2).expanduser().resolve(strict=True)

    if dir1.samefile(dir2):
        raise OSError(f"compare_all inputs are the same directory: {dir1}")

    rtol = 1e-5
    rtolN = rtol
    rtolT = rtol
    rtolJ = rtol
    rtolV = rtol
    atol = 1e-8
    atolN = 1e9
    atolT = 100
    atolJ = 1e-7
    atolV = 50

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

        if not np.allclose(A["ne"][1], B["ne"][1], rtolN, atolN, True):
            errs += 1
            logging.error(f"Ne {st}   {abs(A['ne'][1] - B['ne'][1]).max().item():.3e}")

        for k in ("v1", "v2", "v3"):
            if not np.allclose(A[k][1], B[k][1], rtolV, atolV, True):
                errs += 1
                logging.error(f"{k} {st}   {abs(A[k][1] - B[k][1]).max().item():.3e}")

        for k in ("Ti", "Te"):
            if not np.allclose(A[k][1], B[k][1], rtolT, atolT, True):
                errs += 1
                logging.error(f"{k} {st}   {abs(A[k][1] - B[k][1]).max().item():.3e}")

        for k in ("J1", "J2", "J3"):
            if not np.allclose(A[k][1], B[k][1], rtolJ, atolJ, True):
                errs += 1
                logging.error(f"{k} {st}  {abs(A[k][1] - B[k][1]).max().item():.3e}")

        # %% assert time steps have unique output (earth always rotating...)
        if i > 1:
            if np.allclose(ref["ne"][1], A["ne"][1], rtol, atol):
                errs += 1
                logging.error(f"Ne {st} too similar to prior step")

            for k in ("v1", "v2", "v3"):
                if np.allclose(ref[k][1], A[k][1], rtol, atol):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        if i == 3:
            for k in ("Ti", "Te"):
                if np.allclose(ref[k][1], A[k][1], rtol, atol):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        if i == 2:
            for k in ("J1", "J2", "J3"):
                if np.allclose(ref[k][1], A[k][1], rtol, atol):
                    errs += 1
                    logging.error(f"{k} {st} too similar to prior step")

        ref.update(A)

    return errs
