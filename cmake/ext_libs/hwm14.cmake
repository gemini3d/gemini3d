# this enables CMake imported target HWM14::HWM14
include(FetchContent)

find_package(hwm14 CONFIG)

if(NOT hwm14_FOUND)
  FetchContent_Declare(HWM14
    GIT_REPOSITORY ${hwm14_git}
    GIT_TAG ${hwm14_tag})

  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
    FetchContent_MakeAvailable(HWM14)
  elseif(NOT hwm14_POPULATED)
    FetchContent_Populate(HWM14)
    add_subdirectory(${hwm14_SOURCE_DIR} ${hwm14_BINARY_DIR})
  endif()
endif()

add_library(HWM14::HWM14 INTERFACE IMPORTED)
set_target_properties(HWM14::HWM14 PROPERTIES
  INTERFACE_LINK_LIBRARIES hwm14)

# to avoid modifying HWM14 source code
foreach(f hwm123114.bin dwm07b104i.dat gd2qd.dat)
  set(p ${hwm14_SOURCE_DIR}/src/hwm14/${f})
  file(COPY ${p} DESTINATION ${PROJECT_BINARY_DIR})
  install(FILES ${p} DESTINATION bin)
endforeach()
