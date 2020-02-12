#!/usr/bin/env python
"""
runs a job
"""
from pathlib import Path
import argparse
import gemini.job


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('outdir', help='simulation output directory')
    p.add_argument('-mpiexec', help='path to desired mpiexec executable')
    p.add_argument('-gemexe', help='path to desired gemini.bin')
    P = p.parse_args()

    config_file = Path(__file__).parent / 'config.nml'

    gemini.job.runner(P.mpiexec, P.gemexe, config_file, P.outdir)
