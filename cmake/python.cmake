# Qualification uses hashlib.file_digest, introduced in Python 3.11.
find_package(Python 3.11 COMPONENTS Interpreter QUIET)
if(NOT Python_FOUND OR Python_VERSION VERSION_LESS 3.11)
  set(Python_Interpreter_FOUND FALSE)
endif()

# Refresh on every configure, including a changed or removed interpreter.
foreach(module IN ITEMS numpy h5py scipy)
  string(TOUPPER "${module}" upper)
  set(${upper}_FOUND FALSE CACHE BOOL "Python ${module} importable" FORCE)
  if(Python_Interpreter_FOUND)
    execute_process(COMMAND "${Python_EXECUTABLE}" -c "import ${module}; print(${module}.__version__)"
      RESULT_VARIABLE status OUTPUT_VARIABLE version ERROR_VARIABLE error
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(status STREQUAL "0")
      set(${upper}_FOUND TRUE CACHE BOOL "Python ${module} importable" FORCE)
      message(STATUS "Python ${module}: ${version}")
    else()
      message(VERBOSE "Python ${module} import failed: ${error}")
    endif()
  endif()
endforeach()

set(missing_python_dependencies)
if(NOT Python_Interpreter_FOUND)
  list(APPEND missing_python_dependencies "Python >= 3.11 interpreter")
endif()
foreach(module IN ITEMS numpy h5py scipy)
  string(TOUPPER "${module}" upper)
  if(NOT ${upper}_FOUND)
    list(APPEND missing_python_dependencies "${module}")
  endif()
endforeach()
if(missing_python_dependencies)
  if(gemini3d_require_qualification)
    message(FATAL_ERROR "Required qualification dependencies unavailable: ${missing_python_dependencies}")
  endif()
  message(STATUS "Optional qualification dependencies unavailable: ${missing_python_dependencies}; dependent tests are not registered")
endif()

if(gemini3d_python)
  if(NOT Python_Interpreter_FOUND OR NOT NUMPY_FOUND OR NOT H5PY_FOUND)
    message(FATAL_ERROR "gemini3d_python requires Python >= 3.11, numpy and h5py")
  endif()
  execute_process(COMMAND "${Python_EXECUTABLE}" -c "import gemini3d; print(gemini3d.__version__)"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE status OUTPUT_VARIABLE version ERROR_VARIABLE error
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT status STREQUAL "0")
    message(FATAL_ERROR "Failed to get PyGemini version: ${error}")
  endif()
endif()
