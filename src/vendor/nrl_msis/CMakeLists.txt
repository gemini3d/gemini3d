# We assume MSISE00 is always available, but MSIS 2.0 might not be available.

# ensure no interaction with previous settings
set(CMAKE_REQUIRED_LIBRARIES)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_FLAGS)

include(CheckFortranCompilerFlag)

# --- MSISE00 and MSIS 2.0 require legacy workaround due to non-standard Fortran code
check_fortran_compiler_flag(-w nowarn)
set(msis_flags)
if(nowarn)
  set(msis_flags -w)
endif()

# "static" to help avoid missing runtime library issues when used from Matlab or Python

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(static_flag -static)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel AND NOT WIN32)
  # this test doesn't fail on Windows, even though the flag is not for Windows
  set(static_flag -static-intel)
else()
  set(static_flag)
endif()
if(static_flag)
  check_fortran_compiler_flag(${static_flag} static_ok)
endif()
if(NOT static_ok)
  set(static_flag)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Gfortran >= 8 need -Wno-pedantic to allow mismatched array size inhernet to MSIS.
  # "-w" doesn't disable pedantic
  list(APPEND msis_flags -std=legacy -Wno-pedantic -fno-implicit-none -Wno-error=array-bounds -fcheck=no-all)
endif()

# --- MSISE00 is always available

add_library(msis00mod msis00_gfortran.f msise00_data.f)
# keep msis00 NOT an OBJECT to avoid having to link it explicitly everywhere
target_compile_options(msis00mod PRIVATE ${msis_flags})
target_include_directories(msis00mod INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(msis00mod PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)


# --- build abstract MSIS interface

add_library(msis_ifc msis_interface.f90)
# msis_ifc not OBJECT to avoid linking issues on CMake 3.13
target_link_libraries(msis_ifc PRIVATE msis00mod)
target_include_directories(msis_ifc INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(msis_ifc PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# --- optional MSIS 2.0

if(msis20)
  include(${PROJECT_SOURCE_DIR}/cmake/msis2.cmake)
  target_link_libraries(msis_ifc PRIVATE msis2)
  target_compile_options(msis2 PRIVATE ${msis_flags})
else(msis20)
  target_sources(msis_ifc PRIVATE msis2_dummy.f90)
endif(msis20)

# --- for setting up an equilibrium simulation

add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis_ifc const h5fortran::h5fortran)
target_compile_options(msis_setup PRIVATE ${static_flag})
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
