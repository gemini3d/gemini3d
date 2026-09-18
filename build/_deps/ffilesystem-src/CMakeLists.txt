cmake_minimum_required(VERSION 3.25...4.3)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "Please use out of source build like:
  cmake -Bbuild")
endif()

project(ffilesystem
LANGUAGES C CXX
VERSION 6.4.2
)

include(GNUInstallDirs)
include(CheckSourceCompiles)
include(CheckIncludeFile)

set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
# This project isn't using C++ modules; avoid bugs in module scanning CMake 3.28, 3.29 with GCC >= 14.

include(options.cmake)

if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_IGNORE_PREFIX_PATH $ENV{CONDA_PREFIX})
  list(APPEND CMAKE_IGNORE_PATH $ENV{CONDA_PREFIX}/bin $ENV{CONDA_PREFIX}/lib $ENV{CONDA_PREFIX}/Library $ENV{CONDA_PREFIX}/include)
  # need CMAKE_IGNORE_PATH to ensure env var PATH doesn't interfere
  # despite CMAKE_IGNORE_PREFIX_PATH
endif()

if(ffilesystem_fortran)
  enable_language(Fortran)

  if(ffilesystem_iso_fortran_binding)
    include(cmake/iso_fortran_binding.cmake)
  endif()
endif()



# --- RPATH handling for transitive dependencies
if(ffilesystem_ENABLE_RPATH AND BUILD_SHARED_LIBS AND NOT WIN32)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  set(CMAKE_BUILD_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  #set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  # BUILD_WITH_INSTALL_RPATH makes exe build error can't find libffilesystem.*
endif()

if(ffilesystem_fortran AND NOT DEFINED HAVE_F03TYPE)
  if(ffilesystem_trace)
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/cmake/f03type.f90 _src)
    check_source_compiles(Fortran "${_src}" HAVE_F03TYPE)
  else()
    set(HAVE_F03TYPE true)
  endif()
endif()

# C++20 enables nice but optional debugging / UX
# https://cmake.org/cmake/help/latest/prop_tgt/CXX_STANDARD.html
# must at least set C++17 or build can fail from parent projects
# as some compilers default to less than C++17.
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 17)
  if("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 20)
  else()
    set(CMAKE_CXX_STANDARD 17)
  endif()
endif()

message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION} C++${CMAKE_CXX_STANDARD} CMake ${CMAKE_VERSION} ${CMAKE_SYSTEM_NAME} Toolchain ${CMAKE_TOOLCHAIN_FILE}")
if(ANDROID)
  message(STATUS "Android ABI: ${ANDROID_ABI} Level: ${ANDROID_PLATFORM} NDK: ${ANDROID_NDK} CMAKE_ANDROID_STL_TYPE: ${CMAKE_ANDROID_STL_TYPE}")
endif()


message(VERBOSE
"CMAKE_AR ${CMAKE_AR}
CMAKE_RANLIB ${CMAKE_RANLIB}
CMAKE_STRIP ${CMAKE_STRIP}"
)
message(VERBOSE
"CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}
CMAKE_CXX_COMPILER_LINKER_ID: ${CMAKE_CXX_COMPILER_LINKER_ID}
CMAKE_CXX_COMPILER_LINKER_ID: ${CMAKE_CXX_COMPILER_LINKER_ID}
CMAKE_Fortran_COMPILER_LINKER_ID: ${CMAKE_Fortran_COMPILER_LINKER_ID}"
)

include(cmake/stdcpp.cmake)

if(ffilesystem_stdcpp_version)
  message(STATUS "C++ standard library: ${ffilesystem_stdcpp_version}")
endif()

# CPPcheck
if(ffilesystem_cppcheck)
  find_program(cppcheck_exe NAMES cppcheck REQUIRED)
  set(cppcheck_opts --enable=all --inline-suppr --quiet --suppressions-list=${PROJECT_SOURCE_DIR}/cppcheck.supp)
  set(CMAKE_C_CPPCHECK ${cppcheck_exe} --std=c++${CMAKE_CXX_STANDARD} ${cppcheck_opts})
  set(CMAKE_CXX_CPPCHECK ${cppcheck_exe} --std=c++${CMAKE_CXX_STANDARD} ${cppcheck_opts})
endif()

# --- filesystem library

set(_src
src/c.cpp
src/absolute.cpp
src/copy.cpp
src/extra/compiler.cpp
src/extra/cygwin.cpp
src/disk.cpp
src/equivalent.cpp
src/env.cpp
src/executable.cpp
src/home.cpp
src/inquire.cpp
src/lang.cpp
src/lexical.cpp
src/limits.cpp
src/extra/locale.cpp
src/log.cpp
src/memory.cpp
src/mkdir.cpp
src/move.cpp
src/normalize.cpp
src/os.c
src/extra/owner.cpp
src/parent.cpp
src/partition.cpp
src/permissions.cpp
src/platform.cpp
src/pure.cpp
src/relative.cpp
src/resolve.cpp
src/size.cpp
src/space.cpp
src/symlink.cpp
src/tempdir.cpp
src/touch.cpp
src/time.cpp
src/ulimit.cpp
src/which.cpp
src/windows.cpp
)

if(ffilesystem_extra)
  list(APPEND _src
src/extra/case.cpp
src/extra/component.cpp
src/extra/exepath.cpp
src/extra/libpath.cpp
src/extra/random.cpp
src/extra/removable.cpp
src/extra/shell.cpp
src/extra/sysctl.cpp
src/extra/uid.cpp
src/extra/uname.cpp
src/extra/winsock.cpp
  )
endif()


add_library(ffilesystem ${_src})

include(cmake/compilers.cmake)

if(ffilesystem_pic)
  include(CheckPIESupported)
  check_pie_supported()
  set(CMAKE_POSITION_INDEPENDENT_CODE true)
endif()

# use separate object for Fortran to work with Visual Studio generator
if(ffilesystem_fortran)

  add_library(ffilesystem_Fortran OBJECT src/fortran/filesystem.F90)

  target_compile_definitions(ffilesystem_Fortran PRIVATE
  $<$<NOT:$<BOOL:${HAVE_F03TYPE}>>:NO_F03TYPE>
  $<$<BOOL:${ffilesystem_HAVE_ISO_FORTRAN_BINDING}>:HAVE_ISO_FORTRAN_BINDING>
  )
  set_property(TARGET ffilesystem_Fortran PROPERTY
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  )

  target_include_directories(ffilesystem PRIVATE
    $<$<BOOL:${ffilesystem_ISO_Fortran_binding_path}>:${ffilesystem_ISO_Fortran_binding_path}>)

  target_sources(ffilesystem PRIVATE $<TARGET_OBJECTS:ffilesystem_Fortran>)

endif()

target_include_directories(ffilesystem PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_compile_definitions(ffilesystem PUBLIC
"$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<BOOL:${HAVE_CXX_FILESYSTEM}>>:HAVE_CXX_FILESYSTEM>"
"$<$<BOOL:${ffilesystem_HAVE_ISO_FORTRAN_BINDING}>:HAVE_ISO_FORTRAN_BINDING>"
"$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<BOOL:${ffilesystem_HAVE_CLOCK_CAST}>>:HAVE_CLOCK_CAST>"
"$<$<AND:$<COMPILE_LANGUAGE:C,CXX>,$<BOOL:${ffilesystem_trace}>>:FS_TRACE>"
"$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<BOOL:${ffilesystem_HAVE_DLADDR}>>:ffilesystem_HAVE_DLADDR>"
"$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<BOOL:${ffilesystem_HAVE_COPY_FILE_RANGE}>>:ffilesystem_HAVE_COPY_FILE_RANGE>"
"$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<BOOL:${ffilesystem_HAVE_STATX}>>:HAVE_STATX>"
"$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<BOOL:${ffilesystem_extra}>>:ffilesystem_extra>"
"$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<BOOL:${ffilesystem_HAVE_GetFileInformationByName}>>:HAVE_GETFILEINFORMATIONBYNAME>"
)
if(WIN32 OR CYGWIN)
  target_compile_definitions(ffilesystem PRIVATE $<$<COMPILE_LANGUAGE:CXX>:FS_DLL_NAME="$<TARGET_FILE_NAME:ffilesystem>">)
endif()

target_compile_options(ffilesystem PRIVATE ${ffilesystem_cxxflags})

target_compile_features(ffilesystem PRIVATE
$<$<BOOL:${CMAKE_CXX_STANDARD}>:cxx_std_${CMAKE_CXX_STANDARD}>
c_std_99
)

set_target_properties(ffilesystem PROPERTIES
PUBLIC_HEADER include/ffilesystem.h
EXPORT_NAME filesystem
VERSION ${PROJECT_VERSION}
WINDOWS_EXPORT_ALL_SYMBOLS ON
)
if(ffilesystem_warning_as_error)
  set_property(TARGET ffilesystem PROPERTY COMPILE_WARNING_AS_ERROR true)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0)
    set_property(TARGET ffilesystem PROPERTY LINK_WARNING_AS_ERROR true)
  endif()
endif()


file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

target_include_directories(ffilesystem PUBLIC
"$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include;${PROJECT_SOURCE_DIR}/include>"
$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ffilesystem PUBLIC
${GNU_stdfs}
${CMAKE_DL_LIBS}
)
# don't put generator expressions in link_libraries, causes error when consuming by FetchContent
if(WIN32)
  target_link_libraries(ffilesystem PRIVATE Secur32 Shlwapi Userenv Advapi32 ws2_32)
endif()
# Secur32: GetUserNameExW
# Shlwapi: PathIsUNC
# Userenv: GetUserProfileDirectoryW
# Advapi32: OpenProcessToken, GetTokenInformation, LookupAccountSidA
# ws2_32: gethostname

if(ffilesystem_extra AND APPLE)
  # for is_removable
  target_link_libraries(ffilesystem PRIVATE "$<LINK_LIBRARY:FRAMEWORK,DiskArbitration,CoreFoundation>")
endif()

get_property(_ffs_defs TARGET ffilesystem PROPERTY COMPILE_DEFINITIONS)
message(VERBOSE "Ffilesystem COMPILE_DEFINITIONS: ${_ffs_defs}")
get_property(_ffs_flags TARGET ffilesystem PROPERTY COMPILE_OPTIONS)
message(VERBOSE "Ffilesystem COMPILE_OPTIONS: ${_ffs_flags}")
get_property(_ffs_libs TARGET ffilesystem PROPERTY LINK_LIBRARIES)
message(VERBOSE "Ffilesystem LINK_LIBRARIES: ${_ffs_libs}")

if(ffilesystem_fortran)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/filesystem.mod TYPE INCLUDE)
endif()

install(TARGETS ffilesystem EXPORT ${PROJECT_NAME}-targets)

add_library(ffilesystem::filesystem ALIAS ffilesystem)

# --- CLI

if(ffilesystem_cli)
  add_subdirectory(app)
endif()

# --- Tests
if(ffilesystem_BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

include(cmake/install.cmake)

include(FeatureSummary)

add_feature_info(Fortran ffilesystem_fortran "Fortran bindings")
add_feature_info(Cpp HAVE_CXX_FILESYSTEM "C++ stdlib filesystem")
add_feature_info(shared BUILD_SHARED_LIBS "Build shared library: enables lib_path()")
add_feature_info(path_t HAVE_F03TYPE "Fortran path_t object-oriented interface")

feature_summary(WHAT ENABLED_FEATURES DESCRIPTION "${PROJECT_NAME} enabled features:")
