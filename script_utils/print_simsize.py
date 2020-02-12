#!/usr/bin/env python3
import gemini
import argparse


p = argparse.ArgumentParser()
p.add_argument("fn", help="path to simsize.dat")
p = p.parse_args()

print(gemini.get_simsize(p.fn))
