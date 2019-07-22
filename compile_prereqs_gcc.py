#!/usr/bin/env python
"""
Installs prereqs for Gemini program
CMake >= 3.13 needed in general
"""
import subprocess
import shutil
import logging
from pathlib import Path
import pkg_resources
from argparse import ArgumentParser
import typing
import sys

# ========= user parameters ======================
BUILDDIR = 'build'
LOADLIMIT = '-l 4'

# Library parameters
FC = 'gfortran'
CC = 'gcc'
CXX = 'g++'
MPIVERSION = '3.1.3'  # OpenMPI 4 doesn't seem to work with ScalaPack?
MPIFN = f'openmpi-{MPIVERSION}.tar.bz2'
MPIURL = f'https://download.open-mpi.org/release/open-mpi/v3.1/{MPIFN}'
MPISHA1 = 'b3c60e2bdd5a8a8e758fd741f9a5bebb84da5e81'
MPIDIR = f'openmpi-{MPIVERSION}'

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
    raise SystemExit('could not find CMake')
cmake_ver = subprocess.check_output([CMAKE, '--version'], universal_newlines=True).split('\n')[0].split(' ')[2]
if pkg_resources.parse_version(cmake_ver) < pkg_resources.parse_version('3.13'):
    raise SystemExit(f'CMake {cmake_ver} is less than minimum required CMake 3.13')


def openmpi(wipe: bool, dirs: typing.Dict[str, Path]):
    from script_utils.meson_file_download import url_retrieve
    from script_utils.meson_file_extract import extract_tar

    install_lib = dirs['prefix'] / MPIDIR
    source_lib = dirs['workdir'] / MPIDIR

    url_retrieve(MPIURL, MPIFN, ('sha1', MPISHA1))
    extract_tar(MPIFN, source_lib)

    nice = ['nice'] if sys.platform == 'linux' else []
    cmd = nice + ['./configure', f'--prefix={install_lib}', f'CC={CC}', f'CXX={CXX}', f'FC={FC}']

    subprocess.check_call(cmd, cwd=source_lib)

    cmd = nice + ['make', '-C', str(source_lib), '-j', LOADLIMIT, 'install']
    subprocess.check_call(cmd)


def lapack(wipe: bool,  dirs: typing.Dict[str, Path]):
    install_lib = dirs['prefix'] / LAPACKDIR
    source_lib = dirs['workdir'] / LAPACKDIR
    build_lib = source_lib / BUILDDIR

    update(source_lib, LAPACKGIT)

    build_lib.mkdir(exist_ok=True)

    cachefile = build_lib / 'CMakeCache.txt'
    if wipe and cachefile.is_file():
        cachefile.unlink()

    subprocess.check_call([CMAKE,
                           f'-DCMAKE_INSTALL_PREFIX={install_lib}',
                           '-B', str(build_lib), '-S', str(source_lib)])

    subprocess.check_call([CMAKE, '--build', str(build_lib),
                           '--parallel', '--target', 'install'])


def scalapack(wipe: bool, dirs: typing.Dict[str, Path]):
    install_lib = dirs['prefix'] / SCALAPACKDIR
    source_lib = dirs['install'] / SCALAPACKDIR
    build_lib = source_lib / BUILDDIR

    update(source_lib, SCALAPACKGIT)

    cachefile = build_lib / 'CMakeCache.txt'
    if wipe and cachefile.is_file():
        cachefile.unlink()

    subprocess.check_call([CMAKE,
                           f'-DCMAKE_INSTALL_PREFIX={install_lib}',
                           f'-DLAPACK_ROOT={Path(dirs["prefix"]).expanduser() / LAPACKDIR}',
                           '-B', str(build_lib), '-S', str(source_lib)])

    subprocess.check_call([CMAKE, '--build', str(build_lib),
                           '--parallel', '--target', 'install'])


def mumps(wipe: bool, dirs: typing.Dict[str, Path]):
    install_lib = dirs['prefix'] / MUMPSDIR
    source_lib = dirs['workdir'] / MUMPSDIR
    build_lib = source_lib / BUILDDIR

    scalapack_lib = install_lib.parent / SCALAPACKDIR

    update(source_lib, MUMPSGIT)

    cachefile = build_lib / 'CMakeCache.txt'
    if wipe and cachefile.is_file():
        cachefile.unlink()

    subprocess.check_call([CMAKE,
                           f'-DCMAKE_INSTALL_PREFIX={install_lib}',
                           f'-DSCALAPACK_ROOT={scalapack_lib}',
                           '-B', str(build_lib), '-S', str(source_lib)])

    subprocess.check_call([CMAKE, '--build', str(build_lib),
                           '--parallel', '--target', 'install'])


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
    p.add_argument('libs', help='libraries to compile (lapack, scalapack, mumps)', nargs='+')
    p.add_argument('-prefix', help='toplevel path to install libraries under', default='~/lib_gcc')
    p.add_argument('-workdir', help='toplevel path to where you keep code repos', default='~/code')
    p.add_argument('-wipe', help='wipe before completely recompiling libs', action='store_true')
    P = p.parse_args()

    dirs = {'prefix': Path(P.prefix).expanduser().resolve(),
            'workdir': Path(P.workdir).expanduser().resolve(),
            }

    if 'openmpi' in P.libs:
        openmpi(P.wipe, dirs)
    if 'lapack' in P.libs:
        lapack(P.wipe, dirs)
    if 'scalapack' in P.libs:
        scalapack(P.wipe, dirs)
    if 'mumps' in P.libs:
        mumps(P.wipe, dirs)
