# Simple heurisitic script to download and extract a CMake release to the current working directory
# for time/space reasons, we exclude documentation and help files from the extraction
# Usage:
#   cmake -Dversion=<desired_version> -P install_cmake.cmake

cmake_minimum_required(VERSION 3.10)

if(version VERSION_LESS 3.20)
  message(FATAL_ERROR "CMake version must be at least 3.20. Specify command like `cmake -Dversion=4.4.3` with the desired version")
endif()

option(CMAKE_TLS_VERIFY "Enable TLS verification for downloads" ON)

if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
  string(TOLOWER "$ENV{PROCESSOR_ARCHITECTURE}" arch)
  if(arch STREQUAL "amd64")
    set(arch "x86_64")
  endif()
else()
  execute_process(COMMAND uname -m OUTPUT_VARIABLE arch OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

set(arc_type ".tar.gz")
if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
  set(os "linux-${arch}")
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
  set(os "macos-universal")
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
  set(os "windows-${arch}")
  set(arc_type ".zip")
else()
  message(FATAL_ERROR "Unsupported operating system: ${CMAKE_HOST_SYSTEM_NAME}")
endif()

set(archive_name "cmake-${version}-${os}${arc_type}")

set(url "https://github.com/Kitware/CMake/releases/download/v${version}/${archive_name}")

# CMake < 3.17 file(DOWNLOAD) requires an absolute path for the destination file
if(CMAKE_VERSION VERSION_LESS 3.19)
  get_filename_component(archive ${archive_name} ABSOLUTE)
else()
  file(REAL_PATH ${archive_name} archive)
endif()

if(NOT EXISTS ${archive_name})
  file(DOWNLOAD ${url} ${archive} STATUS s TLS_VERIFY ${CMAKE_TLS_VERIFY})
  list(GET s 0 code)
  if(NOT code EQUAL 0)
    list(GET s 1 msg)
    message(FATAL_ERROR "Failed to download ${url}: ${msg}")
  endif()
endif()

if(CMAKE_VERSION VERSION_LESS 4.5)
  execute_process(COMMAND tar -xf ${archive_name}  --exclude=doc/ --exclude=man/ --exclude=Help/)
else()
  file(ARCHIVE_EXTRACT INPUT ${archive_name} PATTERNS_EXCLUDE "doc/" "man/" "Help/")
endif()
