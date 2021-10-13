cmake_minimum_required(VERSION 3.20...3.22)

if(WIN32)
  message(FATAL_ERROR "MPICH does not work on Windows. Use MS-MPI or Intel MPI instead.")
endif()

if(NOT version)
  set(version 3.4.2)
endif()

if(NOT prefix)
  set(prefix "~")
endif()
if(CMAKE_VERSION VERSION_LESS 3.21)
  get_filename_component(prefix ${prefix} ABSOLUTE)
else()
  file(REAL_PATH ${prefix} prefix EXPAND_TILDE)
endif()

set(CMAKE_TLS_VERIFY true)

set(stem mpich-${version})
set(archive ${prefix}/${stem}.tar.gz)

set(url http://www.mpich.org/static/downloads/${version}/${stem}.tar.gz)

set(install_dir ${prefix}/${stem})
set(src_dir ${install_dir}/${stem})

if(EXISTS ${archive})
  file(SIZE ${archive} fsize)
endif()

if(NOT EXISTS ${archive} OR "${fsize}" LESS 1000000)
  message(STATUS "download ${url} to ${archive}")
  file(DOWNLOAD ${url} ${archive} INACTIVITY_TIMEOUT 15)
endif()

if(NOT EXISTS ${src_dir}/configure.ac)
  message(STATUS "extracting ${archive} to ${install_dir}")
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${install_dir})
endif()

# MPICH uses non-standard Fortran syntax
if(DEFINED ENV{FC})
  set(FC $ENV{FC})
endif()

if(NOT FC)
  find_program(FC gfortran)
endif()

if(FC MATCHES gfortran)
  execute_process(COMMAND ${FC} -dumpversion
  OUTPUT_VARIABLE FC_VERSION
  COMMAND_ERROR_IS_FATAL ANY
  TIMEOUT 5
  )
  if(FC_VERSION VERSION_GREATER_EQUAL 10)
    set(FFLAGS -fallow-argument-mismatch)
  endif()
endif()

execute_process(COMMAND ./configure --prefix=${install_dir} --with-device=ch3 FFLAGS=${FFLAGS}
WORKING_DIRECTORY ${src_dir}
COMMAND_ERROR_IS_FATAL ANY
)

include(ProcessorCount)
ProcessorCount(N)

execute_process(COMMAND make -C ${src_dir} -j ${N} install
COMMAND_ERROR_IS_FATAL ANY
)
