#!/usr/bin/env python3
import subprocess
import shutil

matlab = shutil.which("matlab")
if not matlab:
    raise ImportError("Matlab not available")

subprocess.check_call([matlab, "-batch", "exit"], timeout=60)
