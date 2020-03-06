#!/usr/bin/env python3
"""
Compare two directories of simulation outputs

the absolute and relative tolerance account for slight IEEE-754 based differences,
including non-associativity of floating-point arithmetic.
these parameters are a bit arbitrary.

% per MZ Oct 17, 2018:
% Ti,Te=1 K
% ne=1e6 m-3
% vi,v2,v3=1 m/s
% J1,J2,J3 = 1e-9

% MZ wants to change what we consider significant...
% Ti,Te=5 K
% ne=1e7 m-3
% vi,v2,v3=2 m/s
% J1,J2,J3 = 1e-9

"""
from argparse import ArgumentParser
import gemini3d

tol = {
    "rtol": 1e-5,
    "rtolN": 1e-5,
    "rtolT": 1e-5,
    "rtolJ": 1e-5,
    "rtolV": 1e-5,
    "atol": 1e-8,
    "atolN": 1e9,
    "atolT": 100,
    "atolJ": 1e-7,
    "atolV": 50,
}


def main():
    p = ArgumentParser()
    p.add_argument("outdir", help="directory to compare")
    p.add_argument("refdir", help="reference directory")
    p.add_argument("-p", "--plot", help="make plots of differences", action="store_true")
    P = p.parse_args()

    errs = gemini3d.compare_all(P.outdir, P.refdir, tol, P.plot)

    if errs:
        raise SystemExit(f"{errs} compare errors")
    else:
        print(f"OK: Gemini output comparison {P.outdir} {P.refdir}")


if __name__ == "__main__":
    main()
