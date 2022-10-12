find_package(Matlab COMPONENTS MAIN_PROGRAM REQUIRED)

find_path(matgemini_SOURCE_DIR
NAMES setup_gemini3d.m
PATHS ${PROJECT_SOURCE_DIR}/../mat_gemini/
HINTS ${MATGEMINI_ROOT} ENV MATGEMINI ENV MATGEMINI_ROOT
REQUIRED
)

cmake_path(CONVERT "${matgemini_SOURCE_DIR};${matgemini_SOURCE_DIR}/matlab-stdlib/" TO_NATIVE_PATH_LIST MATLABPATH)

if(MATGEMINI_FOUND)
  return()
endif()

execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch "run('${matgemini_SOURCE_DIR}/setup.m'), stdlib.fileio.expanduser('~');"
RESULT_VARIABLE ret
ERROR_VARIABLE err
)

if(NOT ret EQUAL 0)
  message(FATAL_ERROR "MatGemini not available:
  ${ret}
  ${err}"
  )
endif()

set(MATGEMINI_FOUND true CACHE BOOL "MatGemini found")
