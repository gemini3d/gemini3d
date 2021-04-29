include(ExternalProject)

find_package(glow CONFIG)
if(glow_FOUND)
  return()
endif()

cmake_minimum_required(VERSION 3.19...${CMAKE_VERSION})

if(NOT DEFINED GLOW_ROOT)
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(GLOW_ROOT ${PROJECT_BINARY_DIR}/glow)
  else()
    set(GLOW_ROOT ${CMAKE_INSTALL_PREFIX})
  endif()
endif()

set(GLOW_LIBRARIES
${GLOW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}glow${CMAKE_STATIC_LIBRARY_SUFFIX}
${GLOW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}msis00${CMAKE_STATIC_LIBRARY_SUFFIX})
set(GLOW_INCLUDE_DIRS ${GLOW_ROOT}/include)

ExternalProject_Add(GLOW
GIT_REPOSITORY ${glow_git}
GIT_TAG ${glow_tag}
UPDATE_DISCONNECTED ${EP_UPDATE_DISCONNECTED}
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${GLOW_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
BUILD_BYPRODUCTS ${GLOW_LIBRARIES}
INACTIVITY_TIMEOUT 30
CONFIGURE_HANDLED_BY_BUILD true
)

file(MAKE_DIRECTORY ${GLOW_INCLUDE_DIRS})
# avoid generate race condition

add_library(glow::glow INTERFACE IMPORTED)
target_link_libraries(glow::glow INTERFACE ${GLOW_LIBRARIES})
target_include_directories(glow::glow INTERFACE ${GLOW_INCLUDE_DIRS})
target_compile_definitions(glow::glow INTERFACE DATADIR="${GLOW_ROOT}/share/data/glow/")

add_dependencies(glow::glow GLOW)
