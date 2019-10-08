include(FetchContent)

FetchContent_Declare(ncarglow
  GIT_REPOSITORY https://github.com/gemini3d/glow.git
  GIT_TAG 4828736
)

FetchContent_GetProperties(ncarglow)

if(NOT ncarglow_POPULATED)
  FetchContent_Populate(ncarglow)
  # builds under _deps/ncarglow-build/
  add_subdirectory(${ncarglow_SOURCE_DIR} ${ncarglow_BINARY_DIR})
endif()
