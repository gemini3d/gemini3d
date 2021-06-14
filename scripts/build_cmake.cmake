#!/usr/bin/env -S cmake -P

# NOTE: most users should use install_cmake.cmake instead.
#
# this script builds and installs a recent CMake version
#
# cmake -P build_cmake.cmake
#
# will install CMake under the user's home directory.
#
# optionally, specify a specific CMake version like:
#   cmake -Dversion="3.13.5" -P install_cmake.cmake
#
# This script can be used to install CMake >= 2.8.12.2 (e.g. for compatibility tests)
# old CMake versions have broken file(DOWNLOAD)--they just "download" 0-byte files.

cmake_minimum_required(VERSION 3.7...${CMAKE_VERSION})

set(CMAKE_TLS_VERIFY true)

if(NOT prefix)
  get_filename_component(prefix ~ ABSOLUTE)
endif()

if(NOT version)
  file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/CMAKE_VERSION version
   REGEX "^([0-9]+\.[0-9]+\.[0-9]+)" LIMIT_INPUT 16 LENGTH_MAXIMUM 16 LIMIT_COUNT 1)
endif()

if(version VERSION_LESS 2.8.12.2)
  set(version 2.8.12.2)
endif()

set(host https://github.com/Kitware/CMake/releases/download/v${version}/)
set(stem cmake-${version})
set(name ${stem}.tar.gz)

function(checkup exe)

get_filename_component(path ${exe} DIRECTORY)
set(ep $ENV{PATH})
if(NOT ep MATCHES ${path})
  message(STATUS "add to environment variable PATH ${path}")
endif()

endfunction(checkup)

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/${stem})



find_program(cmake NAMES cmake PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(cmake)
  message(STATUS "CMake ${version} already at ${cmake}")

  checkup(${cmake})
  return()
endif()

message(STATUS "installing CMake ${version} to ${path}")

set(archive ${prefix}/${name})

if(NOT CMAKE_VERSION VERSION_LESS 3.14)
  if(EXISTS ${archive})
    file(SIZE ${archive} fsize)
    if(fsize LESS 1000000)
      file(REMOVE ${archive})
    endif()
  endif()
endif()

if(NOT EXISTS ${archive})
  set(url ${host}${name})
  message(STATUS "download ${url}")
  file(DOWNLOAD ${url} ${archive} INACTIVITY_TIMEOUT 15)

  if(NOT CMAKE_VERSION VERSION_LESS 3.14)
    file(SIZE ${archive} fsize)
    if(fsize LESS 1000000)
      message(FATAL_ERROR "failed to download ${url}")
    endif()
  endif()
endif()

if(NOT IS_DIRECTORY ${path})
  message(STATUS "extracting to ${path}")
  if(CMAKE_VERSION VERSION_LESS 3.18)
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${prefix})
  else()
    file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${prefix})
  endif()
endif()

file(MAKE_DIRECTORY ${path}/build)

execute_process(
  COMMAND ${CMAKE_COMMAND} .. -DBUILD_TESTING:BOOL=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_USE_OPENSSL:BOOL=ON -DCMAKE_INSTALL_PREFIX:PATH=${path}
  RESULT_VARIABLE err
  WORKING_DIRECTORY ${path}/build)
if(NOT err EQUAL 0)
  message(FATAL_ERROR "failed to configure CMake")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build ${path}/build --parallel
  RESULT_VARIABLE err)
if(NOT err EQUAL 0)
  message(FATAL_ERROR "failed to build CMake")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build ${path}/build --target install
  RESULT_VARIABLE err)
if(NOT err EQUAL 0)
  message(FATAL_ERROR "failed to install CMake")
endif()

find_program(cmake NAMES cmake PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(NOT cmake)
  message(FATAL_ERROR "failed to install CMake from ${archive}")
endif()

checkup(${cmake})
