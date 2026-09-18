cmake_minimum_required(VERSION 3.20)

project(fs_check LANGUAGES C CXX)

option(ffilesystem_fortran "use Fortran" on)

if(ffilesystem_fortran)
  enable_language(Fortran)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
    message(VERBOSE "Fortran Linker: ${CMAKE_Fortran_COMPILER_LINKER_ID}")
  endif()
endif()

enable_testing()

if(GNU_stdfs)
  message(STATUS "fs_check: applying flags ${GNU_stdfs}")
endif()

# this is a standin for the Ffilesystem backend library in C++
add_library(ffilesystem_standin OBJECT lib_fs.cpp)
target_include_directories(ffilesystem_standin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ffilesystem_standin PUBLIC ${GNU_stdfs})
target_compile_features(ffilesystem_standin PUBLIC cxx_std_17)

if(ffilesystem_fortran)
  # this is a Fortran interface for the Ffilesystem standin
  add_library(fs_fortran OBJECT lib_fs.f90)

  target_sources(ffilesystem_standin PRIVATE $<TARGET_OBJECTS:fs_fortran>)

  # this is a Fortran main program using Ffilesystem standin
  add_executable(fortran_main main.f90)
  target_link_libraries(fortran_main PRIVATE ffilesystem_standin)
  set_property(TARGET fortran_main PROPERTY LINKER_LANGUAGE ${ffilesystem_linker_lang})

  add_test(NAME Fortran_fs COMMAND fortran_main)
endif()

# this is a C main program using Ffilesystem standin
add_executable(c_main main.c)
target_link_libraries(c_main PRIVATE ffilesystem_standin)
add_test(NAME C_fs COMMAND c_main)

# when the -L system flag is used, say oneAPI with a too old GCC,
# a main C++ target may fail to link but Fortran main targets may succeed

# this is a C++ main program using Ffilesystem standin
add_executable(cxx_fs main.cpp)
target_link_libraries(cxx_fs PRIVATE ffilesystem_standin)
add_test(NAME Cxx_fs COMMAND cxx_fs)
