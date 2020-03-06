include(${CMAKE_CURRENT_LIST_DIR}/test_compare.cmake)

function(win32_env TESTNAME)
# Windows with HDF5 needs this interesting workaround of
# copying HDF5 dlls to the CMAKE_BINARY_DIR.
# adding the dll path to PATH should have worked, but didn't.

if(WIN32) # for Windows ifort dll
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
  if(NOT EXISTS ${CMAKE_BINARY_DIR}/${_f})
    file(COPY ${HDF5_ROOT}/bin/${_f} DESTINATION ${CMAKE_BINARY_DIR})
  endif()
  endforeach()

endif()

endfunction(win32_env)


function(setup_gemini_test TESTNAME TIMEOUT)

# --- disable tests as needed
if(glow)
  set(glow_disabled false)
else()
  string(FIND ${TESTNAME} "glow" _loc)
  if(NOT _loc EQUAL -1)
    set(glow_disabled true)
  endif()
endif()

# --- setup test
set(_outdir ${CMAKE_CURRENT_BINARY_DIR}/test${TESTNAME})
set(_nml ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test${TESTNAME}/inputs/config.nml)

if(python_disabled)
  add_test(NAME gemini:${TESTNAME}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:gemini.bin> ${_nml} ${_outdir}
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
else()
  add_test(NAME gemini:${TESTNAME}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/meson_run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_nml} ${_outdir})
endif()

# NOTE: don't use REQUIRED_FILES because it won't let file download if not present.
set_tests_properties(gemini:${TESTNAME} PROPERTIES
  TIMEOUT ${TIMEOUT}
  SKIP_RETURN_CODE 77
  RUN_SERIAL true
  DISABLED ${glow_disabled}
)
win32_env(gemini:${TESTNAME})

compare_gemini_output(${TESTNAME})

endfunction(setup_gemini_test)
