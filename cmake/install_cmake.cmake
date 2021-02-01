# this script is to install a recent CMake version
# this handles the most common cases, but doesn't handle corner cases like 64-bit kernel with 32-bit user space
# CMAKE_HOST_SYSTEM_PROCESSOR, CMAKE_HOST_SYSTEM_NAME don't work in CMake script mode
#
# cmake -P install_cmake.cmake
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

get_filename_component(prefix ${prefix} ABSOLUTE)

set(ver 3.19.4)

message(STATUS "installing CMake ${ver} to ${prefix}")

set(host https://github.com/Kitware/CMake/releases/download/v${ver}/)


if(APPLE)
  message(STATUS "please use Homebrew https://brew.sh to install cmake:  'brew install cmake'
  or use Python  'pip install cmake'")
  return()
elseif(UNIX)
  execute_process(COMMAND uname -m OUTPUT_VARIABLE arch OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(arch STREQUAL x86_64)
    set(stem cmake-${ver}-Linux-x86_64)
  elseif(arch STREQUAL aarch64)
    set(stem cmake-${ver}-Linux-aarch64)
  endif()
  set(name ${stem}.tar.gz)
elseif(WIN32)
  # https://docs.microsoft.com/en-us/windows/win32/winprog64/wow64-implementation-details?redirectedfrom=MSDN#environment-variables
  # CMake doesn't currently have binary downloads for ARM64 or IA64
  set(arch $ENV{PROCESSOR_ARCHITECTURE})
  if(arch STREQUAL AMD64)
    set(stem cmake-${ver}-win64-x64)
  elseif(arch STREQUAL x86)
    set(stem cmake-${ver}-win32-x86)
  endif()
  set(name ${stem}.zip)
endif()

if(NOT stem)
  message(FATAL_ERROR "unknown CPU arch ${arch}")
endif()

set(path ${prefix}/${stem})

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

message(STATUS "extracting to ${path}")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${prefix})
else()
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${prefix})
endif()

find_program(cmake NAMES cmake PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
if(cmake)
  get_filename_component(path ${cmake} DIRECTORY)
  message(STATUS "add to environment variable PATH ${path}")
else()
  message(FATAL_ERROR "failed to install CMake from ${archive}")
endif()
