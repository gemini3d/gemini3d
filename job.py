#!/usr/bin/env python
"""
runs a job
"""
import argparse
import gemini.job
import sys


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("config_file", help="path to config*.nml file")
    p.add_argument("out_dir", help="simulation output directory")
    p.add_argument("-mpiexec", help="path to desired mpiexec executable")
    p.add_argument("-gemexe", help="path to desired gemini.bin")
    P = p.parse_args()

    try:
        gemini.job.runner(P.mpiexec, P.gemexe, P.config_file, P.out_dir)
    except FileNotFoundError:
        print("\nA necessary simulation input file was not found."
              "\nThis can mean that the simulation initialization script wasn't run first.\n", file=sys.stderr)
        raise
