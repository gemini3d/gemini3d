# ensure no interaction with previous settings
set(CMAKE_REQUIRED_LIBRARIES)
set(CMAKE_REQURIED_INCLUDES)
set(CMAKE_REQUIRED_FLAGS)

include(CheckFortranCompilerFlag)

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


add_library(msis OBJECT msis00_gfortran.f)
target_compile_options(msis PRIVATE ${msis_flags})

# --- for setting up an equilibrium simulation --


add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis)
target_compile_options(msis_setup PRIVATE ${static_flag})
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

if(${PROJECT_NAME}_BUILD_TESTING)
set(_msis_in ${CMAKE_CURRENT_BINARY_DIR}/msis_in.txt)
set(_msis_lz 4)

if(NOT EXISTS ${_msis_in})
file(WRITE ${_msis_in} "
15123
12345
100, 100, 4, 5
${_msis_lz}
40, 50, 60, 70
0, 20, 40, 60
100,110,150,200
")
endif()

add_test(NAME unit:MSISsetup COMMAND $<TARGET_FILE:msis_setup> ${_msis_in} "-" ${_msis_lz})
set_tests_properties(unit:MSISsetup PROPERTIES
  PASS_REGULAR_EXPRESSION "^[ ]*100.00.*110.00.*150.00.*200.00.*862.8.")
endif(${PROJECT_NAME}_BUILD_TESTING)
