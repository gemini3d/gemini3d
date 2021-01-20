# this enables CMake imported target HWM14::HWM14
include(FetchContent)

FetchContent_Declare(hwm14proj
GIT_REPOSITORY ${hwm14_git}
GIT_TAG ${hwm14_tag}
GIT_SHALLOW true
)

FetchContent_MakeAvailable(hwm14proj)

add_library(HWM14::HWM14 INTERFACE IMPORTED)
set_target_properties(HWM14::HWM14 PROPERTIES
  INTERFACE_LINK_LIBRARIES hwm14)

# TODO: check that .dat files are in PROJECT_BINARY_DIR
