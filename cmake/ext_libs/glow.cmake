# Leave GLOW as FetchContent as we use target properties to pass up DATADIR to src/ionization
include(FetchContent)

find_package(glow CONFIG)
if(glow_FOUND)
  return()
endif()

FetchContent_Declare(GLOW
  GIT_REPOSITORY ${glow_git}
  GIT_TAG ${glow_tag})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  FetchContent_MakeAvailable(GLOW)
elseif(NOT glow_POPULATED)
  FetchContent_Populate(GLOW)
  add_subdirectory(${glow_SOURCE_DIR} ${glow_BINARY_DIR})
endif()

file(MAKE_DIRECTORY ${glow_BINARY_DIR}/include)
# avoid generate race condition
