#!/usr/bin/env python
"""
3D translucent visualization -- can be more helpful than slices for some visualizations

https://docs.enthought.com/mayavi/mayavi/mlab_case_studies.html
https://docs.enthought.com/mayavi/mayavi/auto/mlab_decorations.html
"""
from mayavi import mlab
import argparse

import gemini3d


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fn", help=".dat filename to load directly")
    P = p.parse_args()

    dat = gemini3d.readdata(P.fn)

    for p in ("ne", "v1", "Ti", "Te", "J1", "J2", "J3", "v2", "v3"):
        mlab.pipeline.volume(mlab.pipeline.scalar_field(dat[p][1]))
        mlab.title(p)
        mlab.colorbar()
        mlab.show()
