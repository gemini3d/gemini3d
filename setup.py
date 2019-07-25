#!/usr/bin/env python
import setuptools
import sys

if sys.version_info < (3, 6):
    raise SystemExit("Python >= 3.6 required")

setuptools.setup()
