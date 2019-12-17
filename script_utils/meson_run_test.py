#!/usr/bin/env python3
"""
run test
"""
import argparse
import typing
import subprocess
from pathlib import Path

from meson_cpu_count import get_mpi_count
from meson_file_download import url_retrieve
from meson_file_extract import extract_zip

R = Path(__file__).resolve().parents[1]

zenodo: typing.Dict[str, typing.Any] = {
    '2d_fang': {
        'url': 'https://zenodo.org/record/3477385/files/zenodo2d_fang.zip?download=1',
        'md5': '57d72fd0005247c8eedf122ac4670ad0',
        'dir': R / 'tests/data/zenodo2d_fang',
        'zip': R / 'tests/data/zenodo2d_fang.zip',
    },
    '2d_glow': {
        'url': 'https://zenodo.org/record/3477385/files/zenodo2d_glow.zip?download=1',
        'md5': '557bc6a91d8bf3464abdc5c8784f3042',
        'dir': R / 'tests/data/zenodo2d_glow',
        'zip': R / 'tests/data/zenodo2d_glow.zip',
    },
    '3d_fang': {
        'url': 'https://zenodo.org/record/3477330/files/zenodo3d_fang.zip?download=1',
        'md5': 'cf73d6eb166369c522da7a371492a1ce',
        'dir': R / 'tests/data/zenodo3d_fang',
        'zip': R / 'tests/data/zenodo3d_fang.zip',
    },
    '3d_glow': {
        'url': 'https://zenodo.org/record/3477330/files/zenodo3d_glow.zip?download=1',
        'md5': '3528946525295cc8271aa41bc262d7f1',
        'dir': R / 'tests/data/zenodo3d_glow',
        'zip': R / 'tests/data/zenodo3d_glow.zip',
    },
}


def run_test(testname: str, mpiexec: str, exe: str, nml: str, outdir: str, mpi_count: int = None):
    z = zenodo[testname]
    url_retrieve(z['url'], z['zip'], ('md5', z['md5']))
    extract_zip(z['zip'], z['dir'])

    if not mpi_count:
        mpi_count = get_mpi_count(z['dir'] / 'inputs/simsize.dat')

    # have to get exe as absolute path
    exe_abs = Path(exe).resolve()
    cmd = [mpiexec, '-np', str(mpi_count), str(exe_abs), nml, outdir]
    print(' '.join(cmd))
    subprocess.check_call(cmd, cwd=R)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('testname')
    p.add_argument('mpiexec')
    p.add_argument('exe')
    p.add_argument('nml')
    p.add_argument('outdir')
    p.add_argument('-np', help='force number of MPI images', type=int)
    P = p.parse_args()

    run_test(P.testname, P.mpiexec, P.exe, P.nml, P.outdir, mpi_count=P.np)
