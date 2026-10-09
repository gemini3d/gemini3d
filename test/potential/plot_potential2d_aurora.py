#!/usr/bin/env python3
"""
Created on Mon Apr 20 18:13:14 2026

@author: zettergm
"""

import argparse
import numpy as np
import h5py

from matplotlib.pyplot import figure,pcolormesh,xlabel,ylabel,title,colorbar


def read_potential(fn) -> dict:
    with h5py.File(fn, "r") as f:
        out = {"lx1": f["/lx1"][()],
               "lx2": f["/lx2"][()],
               "lx3": f["/lx3"][()],
               "x1": f["/x1"][:],
               "x2": f["/x2"][:],
               "x3": f["/x3"][:],
               "Phi": f["/Phi"][:],
               "A": f["/A"][:],
               "Ap": f["/Ap"][:],
               "SigH": f["/SigH"][:],
               "B": f["/B"][:],
               "C": f["/C"][:],
               "srcterm": f["/srcterm"][:]}

    assert out["lx1"] == out["x1"].size
    assert out["lx2"] == out["x2"].size
    assert out["lx3"] == out["x3"].size

    return out


def plot_potential(data):

    x2 = data["x2"]
    x3 = data["x3"]
    Phi = data["Phi"]
    srcterm = data["srcterm"]

    figure(figsize=(6, 6))
    pcolormesh(x2, x3, Phi)
    colorbar()
    ylabel("distance [m]")
    xlabel("distance [m]")
    title("2D potential (numerical)")

    #figure()
    #pcolormesh(x2,x3,data["A"])
    #colorbar()
    #title("Pedersen")

    figure()
    pcolormesh(x2,x3,data["Ap"])
    colorbar()
    title("Ap")

    # figure()
    # pcolormesh(x2,x3,data["SigH"])
    # colorbar()

    # figure()
    # pcolormesh(x2,x3,data["B"])
    # colorbar()

    # figure()
    # pcolormesh(x2,x3,data["C"])
    # colorbar()

    figure()
    pcolormesh(x2,x3,srcterm)   # srcterm is -Jpar
    colorbar()
    title("Jpar")

    ## Apparently hdf5 mangles array axes???
    Ex,Ey = np.gradient(-1*Phi.transpose(),x2,x3)   # WHYYYYY
    SigP=data["A"].transpose()
    SigH=data["SigH"].transpose()
    Jx=SigP*Ex-SigH*Ey
    Jy=SigH*Ex+SigP*Ey
    Jxx,_ = np.gradient(Jx,x2,x3)
    _,Jyy = np.gradient(Jy,x2,x3)
    divJ=Jxx+Jyy
    Jpartest=-divJ
    errterm=Jpartest-srcterm
    print("Max error in Jpar:", np.max(np.abs(errterm)))

    figure()
    pcolormesh(x2,x3,Ex.transpose())
    colorbar()
    title("Ex")

    figure()
    pcolormesh(x2,x3,Ey.transpose())
    colorbar()
    title("Ey")

    figure()
    pcolormesh(x2,x3,Jx.transpose())
    colorbar()
    title("Jx")

    figure()
    pcolormesh(x2,x3,Jy.transpose())
    colorbar()
    title("Jy")

    figure()
    pcolormesh(x2,x3,Jpartest.transpose())
    colorbar()
    title("div J")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fn", help="HDF5 file containing potential data")
    parser.add_argument("--noplot", action="store_true", help="Do not display plots")
    args = parser.parse_args()
    fn = args.fn

    data = read_potential(fn)
    if not args.noplot:
        plot_potential(data)
