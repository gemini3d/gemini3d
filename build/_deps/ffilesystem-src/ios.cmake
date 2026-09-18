# Usage:
# cmake --workflow ios

# Minimal iOS toolchain for building static/library targets only.
# This is intended for cross-compile demo builds, not app packaging/signing.
set(CMAKE_SYSTEM_NAME iOS)

if(NOT DEFINED CMAKE_OSX_SYSROOT)
  set(CMAKE_OSX_SYSROOT iphoneos CACHE STRING "Apple SDK (iphoneos or iphonesimulator)")
endif()

if(NOT DEFINED CMAKE_OSX_ARCHITECTURES)
  if(CMAKE_OSX_SYSROOT STREQUAL "iphonesimulator")
    set(CMAKE_OSX_ARCHITECTURES arm64 CACHE STRING "Target architectures")
  else()
    set(CMAKE_OSX_ARCHITECTURES arm64 CACHE STRING "Target architectures")
  endif()
endif()

if(NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0" CACHE STRING "Minimum iOS deployment target")
endif()

# Avoid executable try-compile checks when cross-compiling.
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
