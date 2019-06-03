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

# all libraries installed under $PREFIX/libraryname
PREFIX = '~/lib_intel'

# where you keep your Git repos
WORKDIR = '~/code'

# Library parameters
MUMPSGIT = 'https://github.com/scivision/mumps'
MUMPSDIR = 'mumps'

BUILDDIR = 'build'


# ========= end of user parameters ================

GITEXE = shutil.which('git')


def main(wipe: bool):
    install_lib = Path(PREFIX).expanduser() / MUMPSDIR
    source_lib = Path(WORKDIR).expanduser() / MUMPSDIR

    update(source_lib)

    subprocess.run(['python', 'build.py', 'intel', '-install', str(install_lib)], cwd=source_lib)


def update(path: Path):
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
        subprocess.check_call([GITEXE, 'clone', MUMPSGIT, str(path)])


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-wipe', help='wipe before completely recompiling libs', action='store_true')
    P = p.parse_args()

    main(P.wipe)
