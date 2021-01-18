# We assume MSISE00 is always available, but MSIS 2.0 might not be available.

# ensure no interaction with previous settings
set(CMAKE_REQUIRED_LIBRARIES)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_FLAGS)

include(CheckFortranCompilerFlag)


function(build_msis00)

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

# because of old MSIS coding practice, disable vendored MSIS compiler warnings
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Gfortran >= 8 need -Wno-pedantic to allow mismatched array size inhernet to MSIS.
  # "-w" doesn't disable pedantic
  list(APPEND msis_flags -std=legacy -Wno-pedantic -fno-implicit-none -Wno-error=array-bounds -fcheck=no-all)
endif()

add_library(msis00 msis00_gfortran.f)
# keep msis00 NOT an OBJECT to avoid having to link it explicitly everywhere
target_compile_options(msis00 PRIVATE ${msis_flags})

# --- for setting up an equilibrium simulation --

add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis_ifc)
target_compile_options(msis_setup PRIVATE ${static_flag})
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

endfunction(build_msis00)

#--------------------------

build_msis00()

add_library(msis_ifc OBJECT msis_interface.f90)
target_link_libraries(msis_ifc PRIVATE msis00)
target_include_directories(msis_ifc INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(msis_ifc PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)


if(msis20)
  # must have msis2.cmake in this directory to allow use from other projects
  # that start from this CMakeLists.txt by add_subdirectory(../gemini3d/src/vendor/msis00)
  include(msis2.cmake)
  target_link_libraries(msis_ifc PRIVATE msis2)
else(msis20)
  target_sources(msis_ifc PRIVATE msis2_dummy.f90)
endif(msis20)


if(${PROJECT_NAME}_BUILD_TESTING)
  set(_msis_in ${CMAKE_CURRENT_BINARY_DIR}/msis_in.txt)
  set(_msis_lz 2)

  if(NOT EXISTS ${_msis_in})
    file(WRITE ${_msis_in} "
    123
    12345
    100, 100, 4, 5
    ${_msis_lz}
    40, 70
    0, 60
    100,200
    ")
  endif()

  add_test(NAME unit:MSIS00:setup
    COMMAND $<TARGET_FILE:msis_setup> ${_msis_in} "-" ${_msis_lz} 0)
  set_tests_properties(unit:MSIS00:setup PROPERTIES
    PASS_REGULAR_EXPRESSION "^[ ]*100.00.*200.00.*862.8.")

  if(msis20)
    add_test(NAME unit:MSIS2.0:setup
      COMMAND $<TARGET_FILE:msis_setup> ${_msis_in} "-" ${_msis_lz} 20
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    set_tests_properties(unit:MSIS2.0:setup PROPERTIES
      PASS_REGULAR_EXPRESSION "^[ ]*100.00.*200.00.*863.5.")
  endif(msis20)
endif(${PROJECT_NAME}_BUILD_TESTING)
