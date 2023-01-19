function(check_pygemini)

# Numpy conflicts are a general source of trouble
execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy"
RESULT_VARIABLE ret
OUTPUT_VARIABLE out
ERROR_VARIABLE err
)

if(NOT ret EQUAL 0)
  message(VERBOSE "Problem with Python Numpy:
  ${ret}
  ${out}
  ${err}"
  )

  return()
endif()

# need h5py for Python tests
execute_process(COMMAND ${Python_EXECUTABLE} -c "import h5py,numpy,sys; print(f'Python {sys.version}  Numpy {numpy.__version__} h5py {h5py.__version__}')"
RESULT_VARIABLE ret
OUTPUT_VARIABLE out
ERROR_VARIABLE err
)

if(ret EQUAL 0)
  set(H5PY_FOUND true CACHE BOOL "Python h5py Found")
else()
  message(VERBOSE "Problem with Python h5py:
  ${ret}
  ${out}
  ${err}"
  )

  return()
endif()


if(NOT python)
  return()
endif()

execute_process(COMMAND ${Python_EXECUTABLE} -c "import gemini3d; print(gemini3d.__version__)"
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}  # help avoid Intel Windows false import error due to hdf5 .dll in build dir
RESULT_VARIABLE ret
OUTPUT_VARIABLE out
ERROR_VARIABLE err
OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT ret EQUAL 0)
  message(FATAL_ERROR "Failed to get PyGemini version:
  ${ret}
  ${out}
  ${err}"
  )
endif()

endfunction()


find_package(Python COMPONENTS Interpreter)
if(Python_FOUND AND NOT DEFINED PYGEMINI_FOUND)
  check_pygemini()
endif()
