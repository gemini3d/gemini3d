include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(FetchContent)

if(CMAKE_CXX_STANDARD LESS 20 AND "cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  set(CMAKE_CXX_STANDARD 20)
endif()


if(LINUX)
  check_include_file("linux/magic.h" ffilesystem_HAVE_LINUX_MAGIC)
endif()

# -- CLI tests

if(NOT CMAKE_CROSSCOMPILING)
  include(${CMAKE_CURRENT_LIST_DIR}/cli_tests.cmake)
endif()


set(ffilesystem_fortran_test_flags)
set(ffilesystem_cpp_test_flags)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 13)
  set(ffilesystem_fortran_test_flags "$<$<COMPILE_LANGUAGE:Fortran>:-Wno-maybe-uninitialized>")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(ffilesystem_cpp_test_flags "$<$<COMPILE_LANGUAGE:CXX>:--diag_suppress=nodiscard_doesnt_apply>")
endif()
message(VERBOSE "${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}: Fortran test flags: ${ffilesystem_fortran_test_flags}")
message(VERBOSE "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}: C++ test flags: ${ffilesystem_cpp_test_flags}")

set(test_dirs core symlink)
if(ffilesystem_extra)
  list(APPEND test_dirs extra libpath)
endif()


if(CMAKE_CXX_STANDARD GREATER_EQUAL 20)
  FetchContent_Populate(boost_ut URL ${boost_ut_url})
  set(CMAKE_REQUIRED_INCLUDES ${boost_ut_SOURCE_DIR}/include)
  check_include_file_cxx(boost/ut.hpp HAVE_BOOST_UT)
else()
  unset(HAVE_BOOST_UT CACHE)
endif()


foreach(t IN LISTS test_dirs)

add_subdirectory(${t})

# --- Windows shared DLLs
if(BUILD_SHARED_LIBS AND (CYGWIN OR WIN32))

  get_property(n DIRECTORY ${t} PROPERTY TESTS)

  set_property(TEST ${n} DIRECTORY ${t}
  PROPERTY ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:$<TARGET_FILE_DIR:ffilesystem>"
  )
endif()

endforeach()
