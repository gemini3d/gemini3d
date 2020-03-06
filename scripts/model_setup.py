#!/usr/bin/env python
"""
frontend to simply setup a simulation without running it
usually job.py would be used instead.
"""
import argparse
from gemini.model_setup import model_setup


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("config_file", help="path to config*.nml file")
    p.add_argument("out_dir", help="simulation output directory")
    P = p.parse_args()

    model_setup(P.config_file, P.out_dir)
