#!/usr/bin/env python
"""
Intel Fortran needs only Mumps.
CMake >= 3.12 required
"""
import subprocess
import shutil
import logging
from pathlib import Path
from argparse import ArgumentParser

# ========= user parameters ======================

# all libraries installed under PREFIX/libraryname
PREFIX = '~/lib_gcc'

# where you keep your Git repos
WORKDIR = '~/code'
BUILDDIR = 'build'

# Library parameters
MPIVERSION = '3.1.3'  # MPI 4 doesn't currently work with ScalaPack
MPIFN = f'openmpi-{MPIVERSION}.tar.bz2'
MPIURL = f'https://download.open-mpi.org/release/open-mpi/v3.1/{MPIFN}'
MPISHA1 = 'b3c60e2bdd5a8a8e758fd741f9a5bebb84da5e81'
MPIPREFIX = f'{PREFIX}/openmpi-{MPIVERSION}'

LAPACKGIT = 'https://github.com/Reference-LAPACK/lapack'
LAPACKDIR = 'lapack'

SCALAPACKGIT = 'https://github.com/scivision/scalapack'
SCALAPACKDIR = 'scalapack'

MUMPSGIT = 'https://github.com/scivision/mumps'
MUMPSDIR = 'mumps'

# ========= end of user parameters ================

GITEXE = shutil.which('git')
CMAKE = shutil.which('cmake')
if not CMAKE:
    raise FileNotFoundError('could not find CMake')


def lapack(wipe: bool):
    install_lib = Path(PREFIX).expanduser() / LAPACKDIR
    source_lib = Path(WORKDIR).expanduser() / LAPACKDIR
    build_lib = source_lib / BUILDDIR

    update(source_lib, LAPACKGIT)

    build_lib.mkdir(exist_ok=True)

    cachefile = build_lib / 'CMakeCache.txt'
    if wipe and cachefile.is_file():
        cachefile.unlink()

    subprocess.run([CMAKE,
                    f'-DCMAKE_INSTALL_PREFIX={install_lib}',
                    '-B', str(build_lib), '-S', str(source_lib)])

    subprocess.run([CMAKE, '--build', str(build_lib), '--parallel', '--target', 'install'])


def scalapack(wipe: bool):
    install_lib = Path(PREFIX).expanduser() / SCALAPACKDIR
    source_lib = Path(WORKDIR).expanduser() / SCALAPACKDIR
    build_lib = source_lib / BUILDDIR

    update(source_lib, SCALAPACKGIT)

    cachefile = build_lib / 'CMakeCache.txt'
    if wipe and cachefile.is_file():
        cachefile.unlink()

    subprocess.run([CMAKE,
                    f'-DCMAKE_INSTALL_PREFIX={install_lib}',
                    f'-DLAPACK_ROOT={Path(PREFIX).expanduser() / LAPACKDIR}',
                    '-B', str(build_lib), '-S', str(source_lib)])

    subprocess.run([CMAKE, '--build', str(build_lib), '--parallel', '--target', 'install'])


def mumps(wipe: bool):
    install_lib = Path(PREFIX).expanduser() / MUMPSDIR
    source_lib = Path(WORKDIR).expanduser() / MUMPSDIR

    update(source_lib, MUMPSGIT)

    subprocess.run(['python', 'build.py', 'intel', '-install', str(install_lib)], cwd=source_lib)


def update(path: Path, repo: str):
    """
    Use Git to update a local repo, or clone it if not already existing.

    we use cwd= instead of "git -C" for very old Git versions that might be on your HPC.
    """
    if not GITEXE:
        logging.warning('Git not available, cannot check for library updates')
        return

    if path.is_dir():
        subprocess.check_call([GITEXE, 'pull'], cwd=str(path))
    else:
        subprocess.check_call([GITEXE, 'clone', repo, str(path)])


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-wipe', help='wipe before completely recompiling libs', action='store_true')
    P = p.parse_args()

    lapack(P.wipe)

    scalapack(P.wipe)

    mumps(P.wipe)
