cmake_minimum_required(VERSION 3.20...3.22)

if(WIN32)
  message(FATAL_ERROR "OpenMPI does not work on Windows. Use MS-MPI or Intel MPI instead.")
endif()

file(READ ${CMAKE_CURRENT_LIST_DIR}/versions.json _j)
string(JSON openmpi_version GET ${_j} openmpi)

if(NOT prefix)
  set(prefix "~")
endif()
if(CMAKE_VERSION VERSION_LESS 3.21)
  get_filename_component(prefix ${prefix} ABSOLUTE)
else()
  file(REAL_PATH ${prefix} prefix EXPAND_TILDE)
endif()

set(CMAKE_TLS_VERIFY true)

string(SUBSTRING ${openmpi_version} 0 3 subver)

set(stem openmpi-${openmpi_version})
set(archive ${prefix}/${stem}.tar.bz2)

set(url https://download.open-mpi.org/release/open-mpi/v${subver}/${stem}.tar.bz2)

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

execute_process(COMMAND ./configure --prefix=${install_dir}
WORKING_DIRECTORY ${src_dir}
COMMAND_ERROR_IS_FATAL ANY
)

include(ProcessorCount)
ProcessorCount(N)

execute_process(COMMAND make -C ${src_dir} -j ${N} install
COMMAND_ERROR_IS_FATAL ANY
)
