# this script is to install hwloc, to identify host CPU paramters for gemini3d.run
#
# cmake -P install_hwloc.cmake
# will install hwloc under the user's home directory.

if(NOT prefix)
  if(WIN32)
    set(prefix $ENV{USERPROFILE})
  else()
    set(prefix $ENV{HOME})
  endif()
endif()

set(version 2.4.0)

set(host https://download.open-mpi.org/release/hwloc/v2.4)

if(APPLE)
  message(FATAL_ERROR "instead use https://brew.sh via 'brew install hwloc'")
elseif(WIN32)
  set(stem hwloc-win64-build-${version})
  set(name ${stem}.zip)
else()
  set(stem hwloc-${version})
  set(name ${stem}.tar.bz2)
endif()

set(url ${host}/${name})

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/${stem})

message(STATUS "installing hwloc ${version} to ${path}")

set(archive ${path}/${name})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  if(EXISTS ${archive})
    file(SIZE ${archive} fsize)
    if(fsize LESS 100000)
      file(REMOVE ${archive})
    endif()
  endif()
endif()

if(NOT EXISTS ${archive})
  message(STATUS "download ${url}")
  file(DOWNLOAD ${url} ${archive} TLS_VERIFY ON)
endif()

message(STATUS "extracting to ${path}")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})
else()
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${path})
endif()

if(WIN32)
  find_program(lstopo NAMES lstopo PATHS ${path}/${stem} PATH_SUFFIXES bin NO_DEFAULT_PATH)
  if(lstopo)
    get_filename_component(pathbin ${lstopo} DIRECTORY)
    message(STATUS "add to environment variable PATH ${pathbin}")
    message(STATUS "add environment variable HWLOC_ROOT ${path}/${stem}")
  else()
    message(FATAL_ERROR "failed to install Ninja from ${archive}")
  endif()

  return()
endif()

message(STATUS "Building HWLOC")
if(NOT EXISTS ${path}/${stem}/Makefile)
  execute_process(COMMAND ./configure --prefix=${path} WORKING_DIRECTORY ${path}/${stem})
endif()
execute_process(COMMAND make -j -C ${path}/${stem})
execute_process(COMMAND make install -C ${path}/${stem} RESULT_VARIABLE err)

if(NOT err EQUAL 0)
  message(FATAL_ERROR "Hwloc failed to build via autotools")
endif()

find_program(lstopo NAMES lstopo PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
  if(lstopo)
    get_filename_component(pathbin ${lstopo} DIRECTORY)
    message(STATUS "add to environment variable PATH ${pathbin}")
    message(STATUS "add environment variable HWLOC_ROOT ${path}")
  else()
    message(FATAL_ERROR "failed to install Ninja from ${archive}")
  endif()
  return()
