include(GNUInstallDirs)

message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION} CMake ${CMAKE_VERSION}  Toolchain ${CMAKE_TOOLCHAIN_FILE}")

include(cmake/cpu_count.cmake)

cmake_host_system_information(RESULT host_ramMB QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT host_cpu QUERY PROCESSOR_DESCRIPTION)
math(EXPR host_ramGB "${host_ramMB} / 1000")
message(STATUS "Gemini3D: ${host_ramGB} GB RAM detected on ${CMAKE_HOST_SYSTEM_NAME} with ${host_cpu}.  Detected ${Ncpu} CPU cores.")
if(host_ramGB LESS 2)
  message(STATUS "Minimum RAM is about 2 GB--some tests or simulations may fail due to small memory (RAM)")
endif()


if(realbits EQUAL 32)
  message(VERBOSE " 32-bit real precision")
  set(arith s)
else()
  message(VERBOSE " 64-bit real precision")
  set(realbits 64)
  set(arith d)
endif()

option(dev "developer mode: extra compile warnings")

option(glow "use NCAR GLOW airglow / aurora model" on)

option(hwm14 "use HWM14 neutral winds model")

option(python "Python-based self-checks")
# Matlab checks take much longer than Python, and Python covers much more
option(matlab "Matlab-based self-checks")

option(CMAKE_TLS_VERIFY "verify TLS certificates when downloading data" on)

# append .debug to debug libraries, because the computation speed penalty is so great
set(CMAKE_DEBUG_POSTFIX .debug)

# to make Gemini3D more usable by external programs, put all Fortran .mod generated module files in a single directory.
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
# to avoid race condition with imported targets consumed by parent project, create this directory
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Necessary for shared library with Visual Studio / Windows oneAPI
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true)

# CMake < 3.21 will error on configure without this
if(CMAKE_VERSION VERSION_LESS 3.21)
  get_property(_not_top DIRECTORY PROPERTY PARENT_DIRECTORY)
  if(NOT _not_top)
    set(${PROJECT_NAME}_IS_TOP_LEVEL true)
  endif()
endif()

option(${PROJECT_NAME}_BUILD_TESTING "build Gemini3D tests" ${${PROJECT_NAME}_IS_TOP_LEVEL})

file(GENERATE OUTPUT .gitignore CONTENT "*")
