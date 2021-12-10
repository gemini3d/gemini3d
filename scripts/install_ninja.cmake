#!/usr/bin/env -S cmake -P

# this script is to install a recent Ninja version
#
# cmake -P install_ninja.cmake
# will install Ninja under the user's home directory.

cmake_minimum_required(VERSION 3.20...3.22)

if(NOT prefix)
  set(prefix "~")
endif()

set(CMAKE_TLS_VERIFY true)

if(NOT version)
  file(READ ${CMAKE_CURRENT_LIST_DIR}/versions.json _j)
  string(JSON version GET ${_j} ninja)
endif()

set(host https://github.com/ninja-build/ninja/releases/download/v${version}/)


function(checkup ninja)

cmake_path(GET ninja PARENT_PATH ninja_path)

set(ep $ENV{PATH})
cmake_path(CONVERT "${ep}" TO_CMAKE_PATH_LIST ep NORMALIZE)

if(NOT ${ninja_path} IN_LIST ep)
  message(STATUS "add to environment variable PATH ${ninja_path}")
endif()

if(NOT DEFINED ENV{CMAKE_GENERATOR})
  message(STATUS "add environment variable CMAKE_GENERATOR Ninja")
endif()

endfunction(checkup)

if(APPLE)
  execute_process(COMMAND uname -m
    OUTPUT_VARIABLE arch
    OUTPUT_STRIP_TRAILING_WHITESPACE
    TIMEOUT 5
    COMMAND_ERROR_IS_FATAL ANY)
  if(arch STREQUAL x86_64)
    set(stem ninja-mac)
  endif()
elseif(UNIX)
  execute_process(COMMAND uname -m
    OUTPUT_VARIABLE arch
    OUTPUT_STRIP_TRAILING_WHITESPACE
    TIMEOUT 5
    COMMAND_ERROR_IS_FATAL ANY)
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

if(CMAKE_VERSION VERSION_LESS 3.21)
  get_filename_component(prefix ${prefix} ABSOLUTE)
else()
  file(REAL_PATH ${prefix} prefix EXPAND_TILDE)
endif()
set(path ${prefix}/ninja-${version})

message(STATUS "installing Ninja ${version} to ${path}")

set(archive ${path}/${name})

set(url ${host}${name})
message(STATUS "download ${url} to ${archive}")
file(DOWNLOAD ${url} ${archive} INACTIVITY_TIMEOUT 15)

message(STATUS "extracting to ${path}")
file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})

find_program(ninja
  NAMES ninja
  PATHS ${path}
  PATH_SUFFIXES bin
  NO_DEFAULT_PATH
  REQUIRED)

checkup(${ninja})
