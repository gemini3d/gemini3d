# Usage:
# cmake -B builda --toolchain android.cmake -Dffilesystem_fortran=off
# cmake --build builda

if(NOT ANDROID_NDK)
  if(DEFINED ENV{ANDROID_NDK_HOME})
    set(ANDROID_NDK "$ENV{ANDROID_NDK_HOME}" CACHE PATH "Android NDK root")
  elseif(DEFINED ENV{ANDROID_SDK_ROOT})
    file(GLOB _ndk_candidates "$ENV{ANDROID_SDK_ROOT}/ndk/*")
    list(SORT _ndk_candidates DESCENDING)
    list(LENGTH _ndk_candidates _ndk_len)
    if(_ndk_len GREATER 0)
      list(GET _ndk_candidates 0 _ndk_latest)
      set(ANDROID_NDK "${_ndk_latest}" CACHE PATH "Android NDK root")
    endif()
  elseif(DEFINED ENV{ANDROID_HOME} AND EXISTS "$ENV{ANDROID_HOME}/ndk-bundle/build/cmake/android.toolchain.cmake")
    set(ANDROID_NDK "$ENV{ANDROID_HOME}/ndk-bundle" CACHE PATH "Android NDK root")
  endif()
endif()

if(NOT ANDROID_NDK OR NOT EXISTS "${ANDROID_NDK}/build/cmake/android.toolchain.cmake")
  message(FATAL_ERROR "Set ANDROID_NDK or ANDROID_NDK_HOME to a valid NDK path")
endif()

# Defaults (override via -DANDROID_ABI=... -DANDROID_PLATFORM=...)
set(ANDROID_ABI "arm64-v8a" CACHE STRING "Android ABI (arm64-v8a, armeabi-v7a, x86_64, x86)")
set(ANDROID_PLATFORM "latest" CACHE STRING "Android API level (https://apilevels.com/)")
# set(ANDROID_STL "c++_shared" CACHE STRING "NDK C++ STL")

# Delegate to official NDK toolchain
include("${ANDROID_NDK}/build/cmake/android.toolchain.cmake")
