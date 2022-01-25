#!/usr/bin/env -S cmake -P

# this script is to install a recent CMake version
# this handles the most common cases, but doesn't handle corner cases like 64-bit kernel with 32-bit user space
# CMAKE_HOST_SYSTEM_PROCESSOR, CMAKE_HOST_SYSTEM_NAME don't work in CMake script mode
#
#   cmake -P install_cmake.cmake
# will install CMake under the user's home directory.
#
# optionally, specify a specific CMake version like:
#   cmake -Dversion="3.13.5" -P install_cmake.cmake
#
# old CMake versions have broken file(DOWNLOAD)--they just "download" 0-byte files.

cmake_minimum_required(VERSION 3.19...3.22)

set(CMAKE_TLS_VERIFY true)

if(NOT prefix)
  get_filename_component(prefix ~ ABSOLUTE)
endif()

file(READ ${CMAKE_CURRENT_LIST_DIR}/versions.json _j)

if(version VERSION_LESS 3.1)
  # URL patterns changed in CMake 3.1
  string(JSON version GET ${_j} cmake latest)
endif()

# only major.minor specified -- default to latest release known.
string(LENGTH ${version} L)
if (L LESS 5)  # 3.x or 3.xx
  string(JSON version GET ${_j} cmake ${version})
endif()

string(JSON host GET ${_j} cmake binary)
set(host ${host}v${version}/)


function(checkup exe)

get_filename_component(path ${exe} DIRECTORY)
set(ep $ENV{PATH})
if(NOT ep MATCHES ${path})
  message(STATUS "add to environment variable PATH ${path}")
endif()

endfunction(checkup)


if(APPLE)

if(version VERSION_LESS 3.19)
  set(stem cmake-${version}-Darwin-x86_64)
else()
  set(stem cmake-${version}-macos-universal)
endif()

set(name ${stem}.tar.gz)

elseif(UNIX)

execute_process(COMMAND uname -m
OUTPUT_VARIABLE arch
OUTPUT_STRIP_TRAILING_WHITESPACE
TIMEOUT 5
)

if(arch STREQUAL x86_64)
  if(version VERSION_LESS 3.20)
    set(stem cmake-${version}-Linux-x86_64)
  else()
    set(stem cmake-${version}-linux-x86_64)
  endif()
elseif(arch STREQUAL aarch64)
  if(version VERSION_LESS 3.20)
    set(stem cmake-${version}-Linux-aarch64)
  else()
    set(stem cmake-${version}-linux-aarch64)
  endif()
endif()

set(name ${stem}.tar.gz)

elseif(WIN32)

# https://docs.microsoft.com/en-us/windows/win32/winprog64/wow64-implementation-details?redirectedfrom=MSDN#environment-variables
# CMake doesn't currently have binary downloads for ARM64 or IA64
set(arch $ENV{PROCESSOR_ARCHITECTURE})

if(arch STREQUAL AMD64)
  if(version VERSION_LESS 3.6)
    set(stem cmake-${version}-win32-x86)
  elseif(version VERSION_LESS 3.20)
    set(stem cmake-${version}-win64-x64)
  else()
    set(stem cmake-${version}-windows-x86_64)
  endif()
elseif(arch STREQUAL x86)
  if(version VERSION_LESS 3.20)
    set(stem cmake-${version}-win32-x86)
  else()
    set(stem cmake-${version}-windows-i386)
  endif()
endif()

set(name ${stem}.zip)

endif()


if(NOT stem)
  message(FATAL_ERROR "unknown CPU arch ${arch}.  Try building CMake from source:
    cmake -P ${CMAKE_CURRENT_LIST_DIR}/build_cmake.cmake
  or use Python:
    pip install cmake")
endif()

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/${stem})

find_program(cmake NAMES cmake PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(cmake)
  message(STATUS "CMake ${version} already at ${cmake}")

  checkup(${cmake})
  return()
endif()

message(STATUS "installing CMake ${version} to ${prefix}")

set(archive ${prefix}/${name})

if(EXISTS ${archive})
  file(SIZE ${archive} fsize)
  if(fsize LESS 1000000)
    file(REMOVE ${archive})
  endif()
endif()

if(NOT EXISTS ${archive})
  set(url ${host}${name})
  message(STATUS "download ${url}")
  file(DOWNLOAD ${url} ${archive} INACTIVITY_TIMEOUT 15)

  file(SIZE ${archive} fsize)
  if(fsize LESS 1000000)
    if(fsize EQUAL 0)
      file(REMOVE ${archive})
    endif()
    message(FATAL_ERROR "failed to download ${url}")
  endif()
endif()

message(STATUS "extracting to ${path}")
if(CMAKE_VERSION VERSION_LESS 3.18)
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${prefix})
else()
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${prefix})
endif()

find_program(cmake NAMES cmake PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(NOT cmake)
  message(FATAL_ERROR "failed to install CMake from ${archive}")
endif()

checkup(${cmake})
