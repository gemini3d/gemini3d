cmake_policy(VERSION 3.11)

include(FetchContent)

FetchContent_Declare(ncarglow
  GIT_REPOSITORY https://github.com/space-physics/ncar-glow.git
  GIT_TAG 6eacabf681
)


FetchContent_GetProperties(ncarglow)

if(NOT ncarglow_POPULATED)
  FetchContent_Populate(ncarglow)
  # builds under _deps/ncarglow-build/
  add_subdirectory(${ncarglow_SOURCE_DIR} ${ncarglow_BINARY_DIR})
endif()
