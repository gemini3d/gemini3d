# this enables CMake imported target HWM14::HWM14
include(FetchContent)

if(CMAKE_VERSION VERSION_LESS 3.19)
  message(FATAL_ERROR "HWM14 requires CMake >= 3.19")
endif()

FetchContent_Declare(HWM14
  GIT_REPOSITORY ${hwm14_git}
  GIT_TAG ${hwm14_tag})
FetchContent_MakeAvailable(HWM14)

add_library(HWM14::HWM14 INTERFACE IMPORTED)
set_target_properties(HWM14::HWM14 PROPERTIES
  INTERFACE_LINK_LIBRARIES hwm14)

foreach(f hwm123114.bin dwm07b104i.dat gd2qd.dat)
  set(p ${hwm14_SOURCE_DIR}/src/hwm14/${f})
  file(COPY ${p} DESTINATION ${PROJECT_BINARY_DIR})
  install(FILES ${p} DESTINATION bin)
endforeach()
