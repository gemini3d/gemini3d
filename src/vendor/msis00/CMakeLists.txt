include(CheckFortranCompilerFlag)
check_fortran_compiler_flag(-w nowarn)
set(msis_flags)
if(nowarn)
  set(msis_flags -w)
endif()

# static to help avoid missing runtime library issues

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_REQUIRED_FLAGS -static)
  check_fortran_source_compiles("end" static_ok SRC_EXT f90)

  if(static_ok)
    set(static_flags -static)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set(CMAKE_REQUIRED_FLAGS -static-intel)
  check_fortran_source_compiles("end" static_ok SRC_EXT f90)
  if(static_ok)
    set(static_flags -static-intel)
  endif()
endif()

# because of old MSIS coding practice, disable vendored MSIS compiler warnings
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Gfortran >= 8 need -Wno-pedantic to allow mismatched array size inhernet to MSIS.
  # "-w" doesn't disable pedantic
  list(APPEND msis_flags -std=legacy -Wno-pedantic -fno-implicit-none -Wno-error=array-bounds -fcheck=no-all)
endif()


add_library(msis OBJECT msis00_gfortran.f)
target_compile_options(msis PRIVATE ${msis_flags})

# --- for setting up an equilibrium simulation --

add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis)
target_compile_options(msis_setup PRIVATE ${static_flags})
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
