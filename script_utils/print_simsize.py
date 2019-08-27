#!/usr/bin/env python3
from meson_cpu_count import get_simsize
import argparse
from pathlib import Path


p = argparse.ArgumentParser()
p.add_argument("fn", help="path to simsize.dat")
p = p.parse_args()

fn = Path(p.fn).expanduser().resolve()

print(get_simsize(fn))
