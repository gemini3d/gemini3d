cmake_minimum_required(VERSION 3.13)
project(MSIS LANGUAGES Fortran)

if(NOT msis_src)
  if(DEFINED ENV{GEMINI_ROOT})
    set(msis_src $ENV{GEMINI_ROOT}/src/vendor/msis00)
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../vendor/msis00/)
    set(msis_src ${CMAKE_CURRENT_SOURCE_DIR}/../vendor/msis00/)
  else()
    message(FATAL_ERROR "Please define msis_src path to gemini")
  endif()
endif()

add_subdirectory(${msis_src} ${PROJECT_BINARY_DIR}/msis)
