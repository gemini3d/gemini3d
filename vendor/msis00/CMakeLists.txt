# because of old coding practice, disable vendored MSIS compiler warnings
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Gfortran >= 8 need -Wno-pedantic to allow mismatched array size inhernet to MSIS.
  # "-w" doesn't disable pedantic
  set(oldflags -std=legacy -Wno-pedantic -w -fno-implicit-none -Wno-error=array-bounds -fcheck=no-all)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(WIN32)
    set(oldflags /nowarn /nocheck)
  else()
    set(oldflags -nowarn -nocheck)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
  set(oldflags -m 4)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
  set(oldflags -w)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  set(oldflags -w)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  set(oldflags -w)
endif()


add_library(msis msis00_gfortran.f)
target_compile_options(msis PRIVATE ${oldflags})


# --- for setting up an equilibrium simulation --

add_executable(msis_setup ${CMAKE_SOURCE_DIR}/setup/MSIS00/call_msis_gfortran.f90)
target_link_libraries(msis_setup PRIVATE msis)
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})