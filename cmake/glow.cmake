# CMake >= 3.11

include(FetchContent)

FetchContent_Declare(ncarglow
  GIT_REPOSITORY https://github.com/scivision/ncar-glow.git
  GIT_TAG 72398f3
)


FetchContent_GetProperties(ncarglow)

if(NOT ncarglow_POPULATED)
  FetchContent_Populate(ncarglow)
  # builds under bin/_deps/ncarglow-build/
  add_subdirectory(${ncarglow_SOURCE_DIR} ${ncarglow_BINARY_DIR})
endif()
