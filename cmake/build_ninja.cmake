# this script is to build and install a recent Ninja version
#
# cmake -P build_cmake.cmake
# will install Ninja under the user's home directory.

if(NOT prefix)
  if(WIN32)
    set(prefix $ENV{USERPROFILE})
  else()
    set(prefix $ENV{HOME})
  endif()
endif()

set(ver 1.10.2)

set(host https://github.com/ninja-build/ninja/archive/)
set(name v${ver}.tar.gz)

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/ninja-${ver})

message(STATUS "installing Ninja ${ver} to ${path}")

find_program(ninja NAMES ninja PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(ninja)
  get_filename_component(path ${ninja} DIRECTORY)
  set(ep $ENV{PATH})
  message(STATUS "Ninja ${ver} already at ${ninja}")
  if(NOT ep MATCHES "${path}")
    message(STATUS "add to environment variable PATH  ${path}")
  endif()
  return()
endif()

set(archive ${path}/${name})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  if(EXISTS ${archive})
    file(SIZE ${archive} fsize)
    if(fsize LESS 10000)
      file(REMOVE ${archive})
    endif()
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
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
    file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})
  else()
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${path})
  endif()
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
if(ninja)
  get_filename_component(path ${ninja} DIRECTORY)
  message(STATUS "add to environment variable PATH ${path}")
else()
  message(FATAL_ERROR "failed to install Ninja from ${archive}")
endif()
