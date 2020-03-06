include(${CMAKE_CURRENT_LIST_DIR}/test_compare.cmake)

function(win32_env TESTNAME)
# Windows with HDF5 needs to be added to PATH for each test
if(WIN32) # for Windows ifort dll
  if(HDF5_ROOT)
    message(DEBUG " using HDF5_ROOT to set ${TESTNAME} Path")
    set_tests_properties(${TESTNAME} PROPERTIES ENVIRONMENT "PATH=${HDF5_ROOT}/bin;$ENV{PATH}")
  elseif(DEFINED ENV{HDF5_ROOT})
    message(DEBUG " using ENV{HDF5_ROOT} to set ${TESTNAME} Path")
    set_tests_properties(${TESTNAME} PROPERTIES ENVIRONMENT "PATH=$ENV{HDF5_ROOT}/bin;$ENV{PATH}")
  endif()
endif()

endfunction(win32_env)


function(setup_gemini_test TESTNAME TIMEOUT)

# --- disable tests as needed
set(_is_disabled python_disabled)

if(NOT glow)
  string(FIND ${TESTNAME} "glow" _loc)
  if(NOT _loc EQUAL -1)
    set(_is_disabled true)
  endif()
endif()

# --- setup test
set(_outdir ${CMAKE_CURRENT_BINARY_DIR}/test${TESTNAME})
set(_nml ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test${TESTNAME}/inputs/config.nml)

add_test(NAME gemini:${TESTNAME}
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/meson_run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_nml} ${_outdir})

# NOTE: don't use REQUIRED_FILES because it won't let file download if not present.
set_tests_properties(gemini:${TESTNAME} PROPERTIES
  TIMEOUT ${TIMEOUT}
  SKIP_RETURN_CODE 77
  RUN_SERIAL true
  DISABLED ${_is_disabled}
)
win32_env(gemini:${TESTNAME})

compare_gemini_output(${TESTNAME})

endfunction(setup_gemini_test)
