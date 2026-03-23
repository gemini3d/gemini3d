# Downloads external library (e.g. MSIS, GLOW, etc.) source code to allow building offline
cmake_minimum_required(VERSION 3.21)

include(FetchContent)

if(NOT DEFINED offline_dir)
  set(offline_dir ${CMAKE_CURRENT_LIST_DIR}/../offline)
endif()

file(MAKE_DIRECTORY "${offline_dir}")
if(NOT IS_DIRECTORY "${offline_dir}")
  message(FATAL_ERROR "Failed to create offline directory: ${offline_dir}")
endif()
if(NOT EXISTS "${offline_dir}/.gitignore")
  file(WRITE "${offline_dir}/.gitignore" "*\n")
endif()

file(READ "${CMAKE_CURRENT_LIST_DIR}/../cmake/libraries.json" json)

foreach(lib IN ITEMS ffilesystem glow h5fortran hwm14 msis)
  string(JSON ${lib}_url GET "${json}" "${lib}")

  set(${lib}_src "${offline_dir}/${lib}-src")

  FetchContent_Populate(${lib}_source
    URL ${${lib}_url}
    SOURCE_DIR ${${lib}_src}
  )
endforeach()
