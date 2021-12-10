#!/usr/bin/env -S cmake -P

# this script is to build and install a recent Ninja version
#
# cmake -P build_cmake.cmake
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

set(host https://github.com/ninja-build/ninja/archive/)
set(name v${version}.zip)

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

set(src_dir ${path}/ninja-${version})

message(STATUS "extracting ${archive} to ${path}")
file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})

file(MAKE_DIRECTORY ${src_dir}/build)

execute_process(
  COMMAND ${CMAKE_COMMAND} -S${src_dir} -B${src_dir}/build -DBUILD_TESTING:BOOL=OFF -DCMAKE_BUILD_TYPE=Release --install-prefix=${path}
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
