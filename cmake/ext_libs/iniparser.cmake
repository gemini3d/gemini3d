include(FetchContent)

FetchContent_Declare(INIPARSER
GIT_REPOSITORY ${iniparser_git}
GIT_TAG ${iniparser_tag}
INACTIVITY_TIMEOUT 15
)

FetchContent_MakeAvailable(INIPARSER)

set(_s ${iniparser_SOURCE_DIR}/src)

add_library(iniparser ${_s}/iniparser.c ${_s}/dictionary.c)
target_include_directories(iniparser INTERFACE ${_s})
