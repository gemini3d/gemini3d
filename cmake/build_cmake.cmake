# this script is to build and install a recent CMake version
#
# cmake -P build_cmake.cmake
# will install CMake under the user's home directory.
#
# this script should work from CMake >= 2.8.12.

if(NOT prefix)
  if(WIN32)
    set(prefix $ENV{USERPROFILE})
  else()
    set(prefix $ENV{HOME})
  endif()
endif()

set(ver 3.19.4)

set(host https://github.com/Kitware/CMake/releases/download/v${ver}/)
set(stem cmake-${ver})
set(name ${stem}.tar.gz)

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/${stem})

message(STATUS "installing CMake ${ver} to ${path}")

find_program(cmake NAMES cmake PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(cmake)
  get_filename_component(path ${cmake} DIRECTORY)
  set(ep $ENV{PATH})
  message(STATUS "CMake ${ver} already at ${cmake}")
  if(NOT ep MATCHES "${path}")
    message(STATUS "add to environment variable PATH  ${path}")
  endif()
  return()
endif()

set(archive ${prefix}/${name})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
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
  file(DOWNLOAD ${url} ${archive} TLS_VERIFY ON SHOW_PROGRESS)
endif()

if(NOT IS_DIRECTORY ${path})
  message(STATUS "extracting to ${path}")
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
    file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${prefix})
  else()
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${prefix})
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
if(cmake)
  get_filename_component(path ${cmake} DIRECTORY)
  message(STATUS "add to environment variable PATH ${path}")
else()
  message(FATAL_ERROR "failed to install CMake from ${archive}")
endif()
