#!/usr/bin/env python3
import gemini3d
import argparse


p = argparse.ArgumentParser()
p.add_argument("fn", help="path to simsize.dat")
p = p.parse_args()

print(gemini3d.get_simsize(p.fn))
