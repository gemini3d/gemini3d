#!/usr/bin/env python
"""
Intel Fortran needs only Mumps.
CMake >= 3.12 required
"""
import subprocess
import shutil
import logging
import os
from pathlib import Path
from argparse import ArgumentParser

# ========= user parameters ======================

# all libraries installed under $PREFIX/libraryname
PREFIX = '~/lib_gemini_intel'

# where you keep your Git repos
WORKDIR = '~/code'

# Library parameters
MUMPSGIT = 'https://github.com/scivision/mumps'
MUMPSDIR = 'mumps'

BUILDDIR = 'build'

# FC = 'mpifort'  # Intel 18
FC = 'mpiifort'  # Intel 19
CC = 'mpiicc'
CXX = 'mpicxx'


# ========= end of user parameters ================

GITEXE = shutil.which('git')
CMAKE = shutil.which('cmake')
if not CMAKE:
    raise FileNotFoundError('could not find CMake')


def main(wipe: bool):
    if not os.environ.get('MKLROOT'):
        raise EnvironmentError('must have set MKLROOT by running compilervars.bat or source compilervars.sh before this script.')

    if os.name == 'nt':
        # unlike the plain compilers, the MPI compiler wrappers in Windows require the .bat
        compilers = {'FC': FC + '.bat', 'CC': CC + '.bat', 'CXX': CXX + '.bat'}
    else:
        compilers = {'FC': FC, 'CC': CC, 'CXX': CXX}

    install_lib = Path(PREFIX).expanduser() / MUMPSDIR
    source_lib = Path(WORKDIR).expanduser() / MUMPSDIR
    build_lib = source_lib / BUILDDIR

    update(source_lib)

    cachefile = build_lib / 'CMakeCache.txt'

    if wipe and cachefile.is_file():
        cachefile.unlink()

    if os.name == 'nt':
        wopts = ['-G', 'MinGW Makefiles', '-DCMAKE_SH=CMAKE_SH-NOTFOUND']
    else:
        wopts = []

    # generate
    cmd = [CMAKE] + wopts + [f'-DCMAKE_INSTALL_PREFIX={install_lib}',
                             '-B', str(build_lib), '-S', str(source_lib)]
    print('\n', ' '.join(cmd), '\n')

    ret = subprocess.run(cmd, env=os.environ.update(compilers))
    if ret.returncode:
        raise SystemExit(ret.returncode)

    # build
    subprocess.check_call([CMAKE, '--build', str(build_lib), '--parallel', '--target', 'install'])


def update(path: Path):
    if not GITEXE:
        logging.warning('Git not available, cannot check for library updates')
        return

    if path.is_dir():
        subprocess.check_call([GITEXE, '-C', str(path), 'pull'])
    else:
        subprocess.check_call([GITEXE, 'clone', MUMPSGIT, str(path)])


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-wipe', help='wipe before completely recompiling libs', action='store_true')
    P = p.parse_args()

    main(P.wipe)
