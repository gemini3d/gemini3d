include(FetchContent)

set(iniparser_cmake_args
-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
)

FetchContent_Declare(INIPARSER
GIT_REPOSITORY ${iniparser_git}
GIT_TAG ${iniparser_tag}
CMAKE_ARGS ${iniparser_cmake_args}
CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
)

FetchContent_MakeAvailable(INIPARSER)

set(_s ${iniparser_SOURCE_DIR}/src)

add_library(iniparser ${_s}/iniparser.c ${_s}/dictionary.c)
target_include_directories(iniparser INTERFACE ${_s})
