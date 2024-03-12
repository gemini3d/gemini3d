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

option(${PROJECT_NAME}_BUILD_TESTING "build Gemini3D tests" ${PROJECT_IS_TOP_LEVEL})

option(CMAKE_TLS_VERIFY "verify TLS certificates when downloading data" on)

# append .debug to debug libraries, because the computation speed penalty is so great
set(CMAKE_DEBUG_POSTFIX .debug)

# to make Gemini3D more usable by external programs, put all Fortran .mod generated module files in a single directory.
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
# to avoid race condition with imported targets consumed by parent project, create this directory
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Necessary for shared library with Visual Studio / Windows oneAPI
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true)

if(PROJECT_IS_TOP_LEVEL AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/local" CACHE PATH "..." FORCE)
endif()

# --- CMAKE_PREFIX_PATH auto-detection
if(NOT DEFINED CMAKE_PREFIX_PATH AND DEFINED ENV{CMAKE_PREFIX_PATH})
  set(CMAKE_PREFIX_PATH "$ENV{CMAKE_PREFIX_PATH}")
endif()

if(DEFINED CMAKE_PREFIX_PATH)
  list(LENGTH CMAKE_PREFIX_PATH N)
  if(N EQUAL 1)
    get_filename_component(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" ABSOLUTE)
    if(NOT IS_DIRECTORY "${CMAKE_PREFIX_PATH}")
      message(STATUS "did not find CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
      set(CMAKE_PREFIX_PATH "")  # to override bad CMAKE_PREFIX_PATH cache or environment variable
    endif()
  endif()
endif()

if(NOT CMAKE_PREFIX_PATH)
  get_filename_component(home "~" ABSOLUTE)
  string(TOLOWER ${CMAKE_Fortran_COMPILER_ID} fid)

  if(IS_DIRECTORY ${home}/libgem_${fid})
    set(CMAKE_PREFIX_PATH ${home}/libgem_${fid} CACHE PATH "prefix path for gemini3d/external libs")
    message(STATUS "Auto-selecting CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
  endif()
endif()

# check that gemini3d/external libraries are installed
set(need_gemext "
Gemini3D requires several external libraries that are one-time installed via this procedure.
'~/libgem' directory is an arbitrary location.

  git clone https://github.com/gemini3d/external
  cmake -S external -B external/build -DCMAKE_INSTALL_PREFIX=~/libgem
  cmake --build external/build

Now, configure and build Gemini3D (from gemini3d/ directory) like this:

  cmake -B build -DCMAKE_PREFIX_PATH=~/libgem
  cmake --build build --parallel
")

file(GENERATE OUTPUT .gitignore CONTENT "*")
