#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# This is for Mac or Linux computer that uses system libraries (apt install, MacPorts, HomeBrew, LinuxBrew, etc.)


MPIPREFIX=
LAPACKPREFIX=
SCALAPACKPREFIX=
MUMPSPREFIX=

OPTS=
#OPTS="-DUSEGLOW=true"


. script_utils/check.sh

. script_utils/build.sh
