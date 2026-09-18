# Do not use find_package(GTest), in general it finds some other library's GTest that may have an incompatible ABI.
include(FetchContent)

set(FETCHCONTENT_QUIET OFF)
set(INSTALL_GTEST OFF)

FetchContent_Declare(GTest
URL "${ffilesystem_googletest_url}"
FIND_PACKAGE_ARGS CONFIG
)

FetchContent_MakeAvailable(GTest)

if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  foreach(t IN ITEMS gtest gtest_main gmock gmock_main)
    if(WIN32)
      # necessary since GoogleTest injects /WX blindly and that fails builds with modern IntelLLVM.
      # https://learn.microsoft.com/en-us/cpp/build/reference/compiler-option-warning-level
      target_compile_options(${t} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/WX->)
    else()
      # necessary to avoid
      # error: unknown warning option '-Wno-implicit-float-size-conversion' [-Werror,-Wunknown-warning-option]
      target_compile_options(${t} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wno-error=unknown-warning-option;-Wno-error=character-conversion>")
    endif()
  endforeach()
endif()

include(GoogleTest)

# When cross-compiling without a target runner (emulator/device), gtest_discover_tests
# executes the test binary on the host, which fails for e.g. Android NDK builds.
# Set CROSSCOMPILING_EMULATOR_AVAILABLE=ON to opt into gtest_discover_tests even when
# CMAKE_CROSSCOMPILING is true (e.g. when an adb wrapper is configured as the runner).
if(CMAKE_CROSSCOMPILING)
  option(CROSSCOMPILING_EMULATOR_AVAILABLE
    "Target runner/emulator is available; use gtest_discover_tests even when cross-compiling" OFF)
else()
  set(CROSSCOMPILING_EMULATOR_AVAILABLE ON)
endif()

# Helper macro: use gtest_discover_tests on native builds or when a runner is available,
# fall back to gtest_add_tests (source-parsing, no host execution) otherwise.
# Skip code 77 is reserved by adb_run.sh for "no device/emulator" cases.
macro(ffs_gtest_register target_name)
  if(CROSSCOMPILING_EMULATOR_AVAILABLE)
    gtest_discover_tests(${target_name} ${ARGN} TEST_LIST _ffs_gtest_list)
  else()
    gtest_add_tests(TARGET ${target_name} TEST_LIST _ffs_gtest_list)
  endif()

  if(_ffs_gtest_list)
    set_tests_properties(${_ffs_gtest_list} PROPERTIES SKIP_RETURN_CODE 77)
  endif()
endmacro()
