#!/usr/bin/env python
"""
Compare two directories of simulation outputs
"""
from gemini.output_compare import compare_all
from argparse import ArgumentParser


def main():
    p = ArgumentParser()
    p.add_argument("dir1")
    p.add_argument("dir2")
    P = p.parse_args()

    compare_all(P.dir1, P.dir2)


if __name__ == "__main__":
    main()
