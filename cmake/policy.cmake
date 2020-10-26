if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "use cmake -B build or similar to avoid building in-source, which is messy")
endif()

set(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo;Debug" CACHE STRING "Build type selections" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "default build type")
endif()

if(CMAKE_VERSION VERSION_EQUAL 3.19.0-rc1)
  message(FATAL_ERROR "CMake 3.19.0-rc1 has breaking bugs for any project. Please use a different CMake version.")
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.20)
  # explicit source file extensions
  cmake_policy(SET CMP0115 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
  # make missing imported targets fail immediately
  cmake_policy(SET CMP0111 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  # saner ALIAS target policies
  cmake_policy(SET CMP0107 NEW)
  cmake_policy(SET CMP0108 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
  cmake_policy(SET CMP0099 NEW)
endif()
