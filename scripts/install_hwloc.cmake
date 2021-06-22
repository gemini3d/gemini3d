# this script is to install hwloc, to identify host CPU parameters for gemini3d.run
#
# cmake -P install_hwloc.cmake
# will install hwloc under the user's home directory.

cmake_minimum_required(VERSION 3.19...${CMAKE_VERSION})

set(CMAKE_TLS_VERIFY true)

if(NOT prefix)
  get_filename_component(prefix ~ ABSOLUTE)
endif()

set(version 2.5.0)

string(SUBSTRING ${version} 0 3 subver)

set(host https://download.open-mpi.org/release/hwloc/v${subver}/)

if(APPLE)
  find_program(brew
    NAMES brew
    PATHS /usr/local /opt/homeebrew
    PATH_SUFFIXES bin)

  if(brew)
    execute_process(COMMAND ${brew} install hwloc
      COMMAND_ERROR_IS_FATAL ANY)
  endif()

  return()
endif()

if(WIN32)
  set(arch $ENV{PROCESSOR_ARCHITECTURE})
  if(arch STREQUAL AMD64)
    set(stem hwloc-win64-build-${version})
    set(sha256 b64f5ebe534d1ad57cdd4b18ab4035389b68802a97464c1295005043075309ea)
    set(name ${stem}.zip)
  else()
    message(FATAL_ERROR "HWloc binaries provided for x86_64 only. May need to build HWloc from source.")
  endif()
else()
  set(stem hwloc-${version})
  set(sha256 a9cf9088be085bdd167c78b73ddf94d968fa73a8ccf62172481ba9342c4f52c8)
  set(name ${stem}.tar.bz2)
endif()

set(url ${host}${name})

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/${stem})

message(STATUS "installing hwloc ${version} to ${path}")

set(archive ${path}/${name})

if(EXISTS ${archive})
  file(SIZE ${archive} fsize)
  if(fsize LESS 100000)
    file(REMOVE ${archive})
  endif()
endif()

if(NOT EXISTS ${archive})
  message(STATUS "download ${url}")
  file(DOWNLOAD ${url} ${archive}
    INACTIVITY_TIMEOUT 15
    EXPECTED_HASH SHA256=${sha256})
endif()


function(check_hwloc)

find_program(lstopo
  NAMES lstopo
  PATHS ${path}
  PATH_SUFFIXES bin
  NO_DEFAULT_PATH
  REQUIRED)

get_filename_component(pathbin ${lstopo} DIRECTORY)
message(STATUS "add to environment variable PATH ${pathbin}")
message(STATUS "add environment variable HWLOC_ROOT ${path}")

endfunction(check_hwloc)


if(WIN32)
  message(STATUS "${archive} => ${path}")
  file(ARCHIVE_EXTRACT INPUT ${archive}
    DESTINATION ${prefix})

  check_hwloc()
  return()
endif()

# --- Non-Windows only

find_program(MAKE_COMMAND NAMES make REQUIRED)

# find tempdir, as cannot extract and install to same directory
# https://systemd.io/TEMPORARY_DIRECTORIES/
if(DEFINED ENV{TMPDIR})
  set(tmpdir $ENV{TMPDIR})
elseif(IS_DIRECTORY /var/tmp)
  set(tmpdir /var/tmp)
elseif(IS_DIRECTORY /tmp)
  set(tmpdir /tmp)
else()
  set(tmpdir ${prefix}/build)
endif()

file(ARCHIVE_EXTRACT INPUT ${archive}
  DESTINATION ${tmpdir})

set(workdir ${tmpdir}/${stem})

message(STATUS "Building HWLOC in ${workdir}")
if(NOT EXISTS ${workdir}/Makefile)
  execute_process(COMMAND ./configure --prefix=${path}
    WORKING_DIRECTORY ${workdir}
    COMMAND_ERROR_IS_FATAL ANY)
endif()
execute_process(COMMAND ${MAKE_COMMAND} -j
  WORKING_DIRECTORY ${workdir}
  COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND ${MAKE_COMMAND} install
  WORKING_DIRECTORY ${workdir}
  COMMAND_ERROR_IS_FATAL ANY)

check_hwloc()
