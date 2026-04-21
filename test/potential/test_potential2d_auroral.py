#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 18:13:14 2026

@author: zettergm
"""

#!/usr/bin/env python3
#from pathlib import Path
#import argparse
import sys
import typing

import numpy as np
import h5py

from matplotlib.pyplot import figure,pcolormesh,xlabel,ylabel,title,colorbar

fn = "/Users/zettergm/Projects/gemini3d/build/test/potential/outdata.h5"
doplot = True

# if not fn.is_file():
#     print(fn, "not found", file=sys.stderr)
#     raise SystemExit(77)

with h5py.File(fn, "r") as f:
    lx1 = f["/lx1"][()]
    lx2 = f["/lx2"][()]
    lx3 = f["/lx3"][()]
    x1 = f["/x1"][:]
    x2 = f["/x2"][:]
    x3 = f["/x3"][:]
    Phi = f["/Phi"][:]
    A = f["/A"][:]
    Ap = f["/Ap"][:]
    SigH = f["/SigH"][:]
    B = f["/B"][:]
    C = f["/C"][:]
    srcterm = f["/srcterm"][:]
assert lx1 == x1.size
assert lx2 == x2.size
assert lx3 == x3.size

if not doplot:
    exit

fg = figure(figsize=(6, 6))
h = pcolormesh(x2, x3, Phi)
colorbar()
ylabel("distance [m]")
xlabel("distance [m]")
title("2D potential (numerical)")

# figure()
# pcolormesh(x2,x3,A)
# colorbar()

# figure()
# pcolormesh(x2,x3,Ap)
# colorbar()

# figure()
# pcolormesh(x2,x3,SigH)
# colorbar()

# figure()
# pcolormesh(x2,x3,B)
# colorbar()

# figure()
# pcolormesh(x2,x3,C)
# colorbar()

figure()
pcolormesh(x2,x3,srcterm)   # srcterm is -Jpar
colorbar()
title("Jpar")

## Apparently hdf5 mangles array axes???
Ex,Ey = np.gradient(-1*Phi.transpose(),x2,x3)   # WHYYYYY
SigP=A.transpose()
SigH=SigH.transpose()
Jx=SigP*Ex-SigH*Ey
Jy=SigH*Ex+SigP*Ey
Jxx,_ = np.gradient(Jx,x2,x3)
_,Jyy = np.gradient(Jy,x2,x3)
divJ=Jxx+Jyy
Jpartest=-divJ
errterm=Jpartest-srcterm

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
