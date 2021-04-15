#!/usr/bin/env -S cmake -P

# this script is to build and install a recent Ninja version
#
# cmake -P build_cmake.cmake
# will install Ninja under the user's home directory.

cmake_minimum_required(VERSION 3.18...${CMAKE_VERSION})

if(NOT prefix)
  get_filename_component(prefix ~ ABSOLUTE)
endif()

set(ver 1.10.2)

set(host https://github.com/ninja-build/ninja/archive/)
set(name v${ver}.tar.gz)

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
endif()

set(src_dir ${path}/ninja-${ver})

if(NOT EXISTS ${src_dir}/ninjaCMakeLists.txt)
  message(STATUS "extracting ${archive} to ${path}")
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})
endif()

file(MAKE_DIRECTORY ${src_dir}/build)

execute_process(
  COMMAND ${CMAKE_COMMAND} .. -DBUILD_TESTING:BOOL=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=${path}
  RESULT_VARIABLE err
  WORKING_DIRECTORY ${src_dir}/build)
if(NOT err EQUAL 0)
  message(FATAL_ERROR "failed to configure Ninja")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build ${src_dir}/build --parallel
  RESULT_VARIABLE err)
if(NOT err EQUAL 0)
  message(FATAL_ERROR "failed to build Ninja")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --install ${src_dir}/build
  RESULT_VARIABLE err)
if(NOT err EQUAL 0)
  message(FATAL_ERROR "failed to install Ninja")
endif()

find_program(ninja NAMES ninja PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(NOT ninja)
  message(FATAL_ERROR "failed to install Ninja from ${archive}")
endif()

checkup(${ninja})
