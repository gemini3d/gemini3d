cmake_minimum_required (VERSION 3.1)
project(MSIS_setup Fortran)

# because of old coding practice, disable vendored MSIS compiler warnings
add_library(msis msis00_gfortran.f)

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  target_compile_options(msis PRIVATE -std=legacy -w -fno-implicit-none)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  target_compile_options(msis PRIVATE -nowarn)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
  target_compile_options(msis PRIVATE -m 4)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
  target_compile_options(msis PRIVATE -w)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  target_compile_options(msis PRIVATE -w)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  target_compile_options(msis PRIVATE -w)
endif()
# --- for setting up an equilibrium simulation --

add_executable(msis_setup ${CMAKE_CURRENT_SOURCE_DIR}/../../setup/MSIS00/call_msis_gfortran.f90)
target_link_libraries(msis_setup PRIVATE msis)
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../)
