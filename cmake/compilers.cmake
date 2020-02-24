include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)

# for CMake >= 3.17 with  cmake -G "Ninja Multi-Config"
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Build type selections" FORCE)
set(CMAKE_NMC_DEFAULT_BUILD_FILE_CONFIG "Release" CACHE STRING "Default Build type" FORCE)

# === compiler setup
# feel free to add more compiler_*.cmake
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_intel.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_gnu.cmake)
endif()

# === OpenMP
# optional, for possible MUMPS speedup
if(openmp)
  find_package(OpenMP COMPONENTS C Fortran)
endif()
# === MPI
# MPI is used throughout Gemini
find_package(MPI REQUIRED COMPONENTS Fortran)
