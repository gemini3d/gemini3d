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

# --- disable tests as needed
if(glow)
  set(gemini_disabled false)
else()
  string(FIND ${TESTNAME} "glow" _loc)
  if(NOT _loc EQUAL -1)
    set(gemini_disabled true)
  endif()
endif()

if(NOT Python3_FOUND)
  set(gemini_disabled true)
endif()

# --- setup test
set(_outdir ${CMAKE_CURRENT_BINARY_DIR}/test${TESTNAME})
set(_nml ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test${TESTNAME}/inputs/config.nml)

# this would need a downloader, but that needs Python to work properly...
# if(python_disabled)
#   add_test(NAME gemini:${TESTNAME}
#       COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:gemini.bin> ${_nml} ${_outdir}
#       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
# else()
add_test(NAME gemini:${TESTNAME}
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_nml} ${_outdir})


# NOTE: don't use REQUIRED_FILES because it won't let file download if not present.
set_tests_properties(gemini:${TESTNAME} PROPERTIES
  TIMEOUT ${TIMEOUT}
  SKIP_RETURN_CODE 77
  RUN_SERIAL true
  DISABLED ${gemini_disabled}
)
win32_env(gemini:${TESTNAME})

compare_gemini_output(${TESTNAME})

endfunction(setup_gemini_test)
