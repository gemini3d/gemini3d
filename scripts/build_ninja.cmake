#!/usr/bin/env -S cmake -P

# this script is to build and install a recent Ninja version
#
# cmake -P build_cmake.cmake
# will install Ninja under the user's home directory.

cmake_minimum_required(VERSION 3.19...3.21)

set(CMAKE_TLS_VERIFY true)

if(NOT prefix)
  get_filename_component(prefix ~ ABSOLUTE)
endif()

if(NOT version)
  file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/NINJA_VERSION version
   REGEX "^([0-9]+\.[0-9]+\.[0-9]+)" LIMIT_INPUT 16 LENGTH_MAXIMUM 16 LIMIT_COUNT 1)
endif()

set(host https://github.com/ninja-build/ninja/archive/)
set(name v${version}.tar.gz)

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
set(path ${prefix}/ninja-${version})

find_program(ninja NAMES ninja PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(ninja)
  message(STATUS "Ninja ${version} already at ${ninja}")
  checkup(${ninja})
  return()
endif()

message(STATUS "installing Ninja ${version} to ${path}")

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
  file(DOWNLOAD ${url} ${archive} INACTIVITY_TIMEOUT 15)
endif()

set(src_dir ${path}/ninja-${version})

if(NOT EXISTS ${src_dir}/ninjaCMakeLists.txt)
  message(STATUS "extracting ${archive} to ${path}")
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})
endif()

file(MAKE_DIRECTORY ${src_dir}/build)

execute_process(
  COMMAND ${CMAKE_COMMAND} .. -DBUILD_TESTING:BOOL=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=${path}
  WORKING_DIRECTORY ${src_dir}/build
  COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND ${CMAKE_COMMAND} --build ${src_dir}/build --parallel
COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND ${CMAKE_COMMAND} --install ${src_dir}/build
COMMAND_ERROR_IS_FATAL ANY)

find_program(ninja
  NAMES ninja
  PATHS ${path}
  PATH_SUFFIXES bin
  NO_DEFAULT_PATH
  REQUIRED)

checkup(${ninja})
