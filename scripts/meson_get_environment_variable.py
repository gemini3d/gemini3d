#!/usr/bin/env python3
"""
Meson can't get environment variables into meson.build,
so we do it with this helper script.

Only one environment variable may be read at a time.

returns nothing '' if variable doesn't exist
"""
import os
import sys

var = os.environ.get(sys.argv[1])
if var:
    print(var, end="")
