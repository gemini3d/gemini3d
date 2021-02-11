include(FetchContent)

if(NOT python)
  unset(PYGEMINI_DIR CACHE)
  return()
endif()

function(check_pygemini)

# Python3::Interpreter does NOT work, use ${Python3_EXECUTABLE}

# __path__ is always iterable: https://docs.python.org/3/reference/import.html#__path__

# Numpy conflicts are a general source of trouble
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy,sys; print(f'Python {sys.version}  Numpy {numpy.__version__}')"
  RESULT_VARIABLE _ok
  TIMEOUT 15)

if(NOT _ok EQUAL 0)
  message(STATUS "Problem with Python Numpy, cannot use PyGemini")
  return()
endif()

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import gemini3d; import gemini3d.run_test; print(gemini3d.__path__[0])"
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}  # help avoid Intel Windows false import error due to hdf5 .dll in build dir
  OUTPUT_VARIABLE PYGEMINI_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _ok
  TIMEOUT 15)

if(_ok EQUAL 0)
  message(STATUS "PyGemini found: ${PYGEMINI_DIR}")
  set(PYGEMINI_DIR ${PYGEMINI_DIR} CACHE PATH "PyGemini path")
endif()

endfunction()

# --- script
find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_FOUND)
  return()
endif()
# keep this in script so it's not scoped in function

FetchContent_Declare(PYGEMINI
  GIT_REPOSITORY ${pygemini_git}
  GIT_TAG ${pygemini_tag})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  FetchContent_MakeAvailable(PYGEMINI)
elseif(NOT pygemini_POPULATED)
  FetchContent_Populate(PYGEMINI)
endif()


if(NOT PYGEMINI_DIR)
  check_pygemini()
endif()

if(NOT PYGEMINI_DIR)

# detect virtualenv
# this is how CMake itself works for FindPython3 in Modules/FindPython/Support.cmake
if(DEFINED ENV{VIRTUAL_ENV}) # OR DEFINED ENV{CONDA_PREFIX})  # intel python HPC admin install
  set(_pip_args)
else()
  set(_pip_args "--user")
endif()

execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install -e ${pygemini_SOURCE_DIR} ${_pip_args}
  RESULT_VARIABLE _ok)

if(_ok EQUAL 0)
  message(STATUS "Setup PyGemini in ${pygemini_SOURCE_DIR}")
else()
  message(STATUS "Problem installing PyGemini with ${Python3_EXECUTABLE}")
  return()
endif()

check_pygemini()

endif(NOT PYGEMINI_DIR)

if(NOT PYGEMINI_DIR)
  message(STATUS "MISSING: PyGemini ${pygemini_git}")
endif()
