#!/usr/bin/env python3
"""
Compare two directories of simulation outputs
"""
import gemini.output_compare
from argparse import ArgumentParser


def main():
    p = ArgumentParser()
    p.add_argument("dir1")
    p.add_argument("dir2")
    P = p.parse_args()

    errs = gemini.output_compare.compare_all(P.dir1, P.dir2)

    if errs:
        raise SystemExit(f"{errs} compare errors")
    else:
        print(f"OK: Gemini output comparison {P.dir1} {P.dir2}")


if __name__ == "__main__":
    main()
