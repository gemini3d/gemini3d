# installs Gemini3D basic prereqs on Linux and MacOS, and Windows with MSYS2
# use by:
#
#  cmake -P scripts/install_prereq.cmake

cmake_minimum_required(VERSION 3.7...3.20)

if(WIN32)
  message(FATAL_ERROR "Please install Gemini prereqs on Windows via MSYS2 Terminal https://www.msys2.org/")
endif()

function(check_ninja)
  find_program(ninja NAMES ninja)
  execute_process(COMMAND ninja --version
    OUTPUT_STRIP_TRAILING_WHITESPACE
    TIMEOUT 5
    OUTPUT_VARIABLE _ver
    COMMAND_ERROR_IS_FATAL ANY)

  if(_ver VERSION_LESS 1.10.0)
    execute_process(COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/install_ninja.cmake
      COMMAND_ERROR_IS_FATAL ANY)
  endif()
endfunction(check_ninja)

execute_process(COMMAND uname -s OUTPUT_VARIABLE id TIMEOUT 5)

if(id MATCHES "MSYS")
  execute_process(COMMAND pacman -S --needed mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-ninja mingw-w64-x86_64-hwloc mingw-w64-x86_64-msmpi mingw-w64-x86_64-hdf5 mingw-w64-x86_64-lapack mingw-w64-x86_64-scalapack
    COMMAND_ERROR_IS_FATAL ANY)
  return()
endif()

if(APPLE)
  find_program(brew
    NAMES brew
    PATHS /usr/local /opt/homebrew
    PATH_SUFFIXES bin)

  if(NOT brew)
    message(FATAL_ERROR "We generally suggest installing Homebrew package manager https://brew.sh")
  endif()

  execute_process(COMMAND ${brew} install gcc ninja cmake hwloc lapack scalapack openmpi hdf5
    COMMAND_ERROR_IS_FATAL ANY)
  return()
endif()

find_program(apt NAMES apt)
if(apt)
  execute_process(COMMAND apt install --no-install-recommends ninja-build gfortran libhwloc-dev libmumps-dev liblapack-dev libscalapack-mpi-dev libopenmpi-dev openmpi-bin libhdf5-dev
    COMMAND_ERROR_IS_FATAL ANY)

  check_ninja()
  return()
endif()

find_program(yum NAMES yum)
if(yum)
  execute_process(COMMAND yum install epel-release gcc-gfortran hwloc-devel MUMPS-openmpi-devel lapack-devel scalapack-openmpi-devel openmpi-devel hdf5-devel
    COMMAND_ERROR_IS_FATAL ANY)

  check_ninja()
  return()
endif()

find_program(pacman NAMES pacman)
if(pacman)
  execute_process(COMMAND pacman -S --needed gcc-fortran ninja hwloc openmpi hdf5 lapack scalapack mumps
    COMMAND_ERROR_IS_FATAL ANY)

  check_ninja()
  return()
endif()
