include(${CMAKE_CURRENT_LIST_DIR}/test_compare.cmake)

function(win32_env TESTNAME)
# Windows with HDF5 needs this interesting workaround of
# copying HDF5 dlls to the CMAKE_BINARY_DIR.
# adding the dll path to PATH should have worked, but didn't.

if(WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL Intel) # for Windows ifort dll
  if(NOT DEFINED HDF5_ROOT)
    if(DEFINED ENV{HDF5_ROOT})
      set(HDF5_ROOT $ENV{HDF5_ROOT})
    else()
      return()
    endif()
  endif()

  # bizarre behavior if not strip quotes
  string(REPLACE "\"" "" HDF5_ROOT ${HDF5_ROOT})
  # This "should" have worked, but the following workaround is needed instead.
  # set_tests_properties(${TESTNAME} PROPERTIES ENVIRONMENT "PATH=${HDF5_ROOT}/bin;$ENV{PATH}")

  foreach(_f hdf5.dll hdf5_f90cstub.dll hdf5_fortran.dll hdf5_hl.dll hdf5_hl_f90cstub.dll hdf5_hl_fortran.dll)
  if(NOT EXISTS ${PROJECT_BINARY_DIR}/${_f})
    message(VERBOSE " ${HDF5_ROOT}/bin/${_f} => ${PROJECT_BINARY_DIR}")
    file(COPY ${HDF5_ROOT}/bin/${_f} DESTINATION ${PROJECT_BINARY_DIR})
  endif()
  endforeach()

endif()

endfunction(win32_env)


function(setup_gemini_test TESTNAME TIMEOUT)

# --- setup test
set(_outdir ${CMAKE_CURRENT_BINARY_DIR}/test${TESTNAME})
set(_nml ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test${TESTNAME}/inputs/config.nml)

if(python_ok AND hdf5)
  add_test(NAME gemini:hdf5:${TESTNAME}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_nml} ${_outdir})

  # NOTE: don't use REQUIRED_FILES because it won't let file download if not present.
  set_tests_properties(gemini:hdf5:${TESTNAME} PROPERTIES
    TIMEOUT ${TIMEOUT}
    SKIP_RETURN_CODE 77
    RUN_SERIAL true
    FIXTURES_REQUIRED MPIMUMPS)
  win32_env(gemini:${TESTNAME})
endif()

if(python_ok AND netcdf)
  add_test(NAME gemini:netcdf:${TESTNAME}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_nml} ${_outdir} -out_format nc)

  set_tests_properties(gemini:netcdf:${TESTNAME} PROPERTIES
    TIMEOUT ${TIMEOUT}
    SKIP_RETURN_CODE 77
    RUN_SERIAL true
    FIXTURES_REQUIRED MPIMUMPS)
endif()


compare_gemini_output(${TESTNAME})

endfunction(setup_gemini_test)
