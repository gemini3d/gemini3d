# this script is to install a recent Ninja version
# this handles the most common cases, but doesn't handle corner cases like 64-bit kernel with 32-bit user space
# CMAKE_HOST_SYSTEM_PROCESSOR, CMAKE_HOST_SYSTEM_NAME don't work in CMake script mode
#
# cmake -P install_ninja.cmake
# will install Ninja under the user's home directory.

if(NOT prefix)
  if(WIN32)
    set(prefix $ENV{USERPROFILE})
  else()
    set(prefix $ENV{HOME})
  endif()
endif()

set(ver v1.10.2)

set(host https://github.com/ninja-build/ninja/releases/download/${ver}/)

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
  message(FATAL_ERROR "unknown CPU arch ${arch}. Try building Ninja from source: 'cmake -P ${CMAKE_CURRENT_LIST_DIR}/build_ninja.cmake'")
endif()

set(name ${stem}.zip)

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/ninja-${ver})

message(STATUS "installing Ninja ${ver} to ${path}")

find_program(ninja NAMES ninja PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(ninja)
  get_filename_component(path ${ninja} DIRECTORY)
  set(ep $ENV{PATH})
  message(STATUS "Ninja ${ver} already at ${ninja}")
  if(NOT ep MATCHES "${ninja}")
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

message(STATUS "extracting to ${path}")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})
else()
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${path})
endif()

find_program(ninja NAMES ninja PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(ninja)
  get_filename_component(path ${ninja} DIRECTORY)
  message(STATUS "add to environment variable PATH ${path}")
else()
  message(FATAL_ERROR "failed to install Ninja from ${archive}")
endif()
