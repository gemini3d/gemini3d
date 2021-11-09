function(check_pygemini)

# Python3::Interpreter does NOT work, use ${Python3_EXECUTABLE}

# __path__ is always iterable: https://docs.python.org/3/reference/import.html#__path__

# Numpy conflicts are a general source of trouble
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy,sys; print(f'Python {sys.version}  Numpy {numpy.__version__}')"
RESULT_VARIABLE _ok
TIMEOUT 15
)

if(NOT _ok EQUAL 0)
  message(STATUS "Problem with Python Numpy, cannot use PyGemini")
  return()
endif()

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import gemini3d; import gemini3d.run_test; print(gemini3d.__path__[0])"
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}  # help avoid Intel Windows false import error due to hdf5 .dll in build dir
OUTPUT_VARIABLE PYGEMINI_DIR
OUTPUT_STRIP_TRAILING_WHITESPACE
RESULT_VARIABLE _ok
TIMEOUT 15
)

if(_ok EQUAL 0)
  message(STATUS "PyGemini found: ${PYGEMINI_DIR}")
  set(PYGEMINI_DIR ${PYGEMINI_DIR} CACHE PATH "PyGemini path")
endif()

endfunction()


find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND AND NOT PYGEMINI_DIR)
  check_pygemini()
endif()
