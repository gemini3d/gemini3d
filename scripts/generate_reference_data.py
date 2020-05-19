#!/usr/bin/env python
"""
Generate reference data from test scenario directory automatically
"""
from pathlib import Path

import gemini3d.job

params = {
        "config_file": Path(P.config_file).expanduser(),
        "out_dir": Path(P.out_dir).expanduser(),
        "mpiexec": P.mpiexec,
        "gemexe": P.gemexe,
        "force": P.force,
        "out_format": P.out_format,
        "cpu_count": P.cpu,
    }

ret = gemini3d.job.runner(params)