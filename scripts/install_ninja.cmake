#!/usr/bin/env -S cmake -P

# this script is to install a recent Ninja version
#
# cmake -P install_ninja.cmake
# will install Ninja under the user's home directory.

cmake_minimum_required(VERSION 3.18...${CMAKE_VERSION})

if(NOT prefix)
  get_filename_component(prefix ~ ABSOLUTE)
endif()

set(ver v1.10.2)

set(host https://github.com/ninja-build/ninja/releases/download/${ver}/)


function(checkup ninja)

get_filename_component(path ${ninja} DIRECTORY)
set(ep $ENV{PATH})
if(NOT ep MATCHES ${path})
  message(STATUS "add to environment variable PATH ${path}")
endif()

if(NOT DEFINED ENV{CMAKE_GENERATOR})
  message(STATUS "add environment variable CMAKE_GENERATOR Ninja")
endif()

endfunction(checkup)

if(APPLE)
  execute_process(COMMAND uname -m OUTPUT_VARIABLE arch OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(arch STREQUAL x86_64)
    set(stem ninja-mac)
  endif()
elseif(UNIX)
  execute_process(COMMAND uname -m OUTPUT_VARIABLE arch OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(arch STREQUAL x86_64)
    set(stem ninja-linux)
  endif()
elseif(WIN32)
  # https://docs.microsoft.com/en-us/windows/win32/winprog64/wow64-implementation-details?redirectedfrom=MSDN#environment-variables
  set(arch $ENV{PROCESSOR_ARCHITECTURE})
  if(arch STREQUAL AMD64)
    set(stem ninja-win)
  endif()
endif()

if(NOT stem)
  message(FATAL_ERROR "unknown CPU arch ${arch}. Try building Ninja from source:
    cmake -P ${CMAKE_CURRENT_LIST_DIR}/build_ninja.cmake")
endif()

set(name ${stem}.zip)

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/ninja-${ver})

find_program(ninja NAMES ninja PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(ninja)
  message(STATUS "Ninja ${ver} already at ${ninja}")
  checkup(${ninja})
  return()
endif()

message(STATUS "installing Ninja ${ver} to ${path}")

set(archive ${path}/${name})

if(EXISTS ${archive})
  file(SIZE ${archive} fsize)
  if(fsize LESS 10000)
    file(REMOVE ${archive})
  endif()
endif()

if(NOT EXISTS ${archive})
  set(url ${host}${name})
  message(STATUS "download ${url}")
  file(DOWNLOAD ${url} ${archive} TLS_VERIFY ON)
  file(SIZE ${archive} fsize)
  if(fsize LESS 10000)
    message(FATAL_ERROR "failed to download ${url}")
  endif()
endif()

message(STATUS "extracting to ${path}")
file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})

find_program(ninja NAMES ninja PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(NOT ninja)
  message(FATAL_ERROR "failed to install Ninja from ${archive}")
endif()

checkup(${ninja})
