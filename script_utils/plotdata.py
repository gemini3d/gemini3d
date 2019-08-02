#!/usr/bin/env python3
from argparse import ArgumentParser
from matplotlib.pyplot import figure, show
from readdata import loadframe


def plotplasma(dat: dict):
    for k, v in dat.items():
        if k == "t":
            t = v
            continue

        ax = figure().gca()
        ax.pcolormesh(v.squeeze())
        ax.set_title(f"{k}  {t}")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("simdir")
    p = p.parse_args()

    dat = loadframe(p.simdir)

    plotplasma(dat)

    show()
