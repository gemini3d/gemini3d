#!/usr/bin/env python3
import argparse
import gemini3d


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("fn", help="path containing simsize.{h5,nc,dat}")
    p.add_argument("-f", "--force", help="force CPU count", type=int)
    P = p.parse_args()

    mpi_count = gemini3d.get_mpi_count(P.fn, P.force)

    # need end='' or you'll have to .strip() in Meson
    print(mpi_count, end="")
