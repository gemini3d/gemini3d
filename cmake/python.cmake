include(${CMAKE_CURRENT_LIST_DIR}/git.cmake)

function(check_pygemini)

# Python3::Interpreter does NOT work, use ${Python3_EXECUTABLE}

# __path__ is always iterable: https://docs.python.org/3/reference/import.html#__path__

find_package(Python3 COMPONENTS Interpreter)

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import gemini3d; print(gemini3d.__path__[0])"
  OUTPUT_VARIABLE PYGEMINI_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _ok
  TIMEOUT 15)

if(_ok EQUAL 0)
  message(STATUS "PyGemini: ${PYGEMINI_DIR}")
  set(python_ok true CACHE BOOL "PyGemini OK")
endif()

endfunction()

# --- script
set(pygemini_dir "${PROJECT_SOURCE_DIR}/../pygemini/")
set(pygemini_url "https://github.com/gemini3d/pygemini")

if(NOT python_ok)
  check_pygemini()
endif()

if(NOT python_ok)
  clone_if_missing(${pygemini_dir} ${pygemini_url})

  execute_process(COMMAND ${Python3_EXECUTABLE} setup.py --user develop
    WORKING_DIRECTORY ${pygemini_dir}
    RESULT_VARIABLE _ok)
  if(NOT _ok EQUAL 0)
    message(WARNING "Problem installing PyGemini. Perhaps try doing this manually: ${pygemini_url}")
    return()
  endif()

  check_pygemini()
endif()

if(NOT python_ok)
  message(WARNING "MISSING: PyGemini ${pygemini_url} many self-tests will not work.")
endif()
