function(check_pygemini)

# Python::Interpreter does NOT work, use ${Python_EXECUTABLE}

# Numpy conflicts are a general source of trouble
execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy,sys; print(f'Python {sys.version}  Numpy {numpy.__version__}')"
RESULT_VARIABLE ret
OUTPUT_VARIABLE out
ERROR_VARIABLE err
TIMEOUT 15
)

if(NOT ret EQUAL 0)
  message(FATAL_ERROR "Problem with Python Numpy, cannot use PyGemini
  ${ret}
  ${out}
  ${err}"
  )
endif()

execute_process(COMMAND ${Python_EXECUTABLE} -c "import gemini3d; print(gemini3d.__version__)"
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}  # help avoid Intel Windows false import error due to hdf5 .dll in build dir
RESULT_VARIABLE ret
OUTPUT_VARIABLE out
ERROR_VARIABLE err
TIMEOUT 15
OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT ret EQUAL 0)
  message(FATAL_ERROR "Failed to get PyGemini version:
  ${ret}
  ${out}
  ${err}"
  )
endif()

set(PYGEMINI_FOUND true CACHE BOOL "PyGemini Found")
set(PYGEMINI_VERSION ${out} CACHE STRING "PyGemini version")

endfunction()


find_package(Python COMPONENTS Interpreter)
if(Python_FOUND AND NOT PYGEMINI_FOUND)
  check_pygemini()
endif()
