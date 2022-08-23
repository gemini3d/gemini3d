# ensure we have a directory for gemini3d/external that is reachable
if(NOT DEFINED CMAKE_PREFIX_PATH AND DEFINED ENV{CMAKE_PREFIX_PATH})
  set(CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
endif()
if(CMAKE_PREFIX_PATH)
  get_filename_component(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ABSOLUTE)
endif()

set(need_gemext "Must specify location of Gemini3D/external libraries like:
cmake -B build -DCMAKE_PREFIX_PATH=~/libgem

where ~/libgem is an arbitrary location you're previously installed libraries to:
https://github.com/gemini3d/external
")

if(NOT CMAKE_PREFIX_PATH)
  message(FATAL_ERROR ${need_gemext})
endif()

# toolchain
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
  set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()
if(CMAKE_TOOLCHAIN_FILE)
  get_filename_component(CMAKE_TOOLCHAIN_FILE ${CMAKE_TOOLCHAIN_FILE} ABSOLUTE)
endif()

if(NOT DEFINED CRAY AND DEFINED ENV{CRAYPE_VERSION})
  set(CRAY true)
endif()

# our cray.cmake from gemini3d/external
if(CRAY)
  if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE ${CMAKE_PREFIX_PATH}/cmake/cray.cmake)
  endif()

  if(NOT EXISTS ${CMAKE_TOOLCHAIN_FILE})
    message(FATAL_ERROR ${need_gemext})
  endif()
endif()
