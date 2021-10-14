include(ProcessorCount)

function(cmake_cpu_count)
  # on ARM e.g. Raspberry Pi, the usually reliable cmake_host_system_info gives 1 instead of true count
  # fallback to less reliable ProcessorCount which does work on Raspberry Pi.
  ProcessorCount(_ncount)
  cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)

  if(Ncpu EQUAL 1 AND _ncount GREATER 0)
    set(Ncpu ${_ncount})
  endif()

  set(Ncpu ${Ncpu} PARENT_SCOPE)

endfunction(cmake_cpu_count)
cmake_cpu_count()

cmake_host_system_information(RESULT host_ramMB QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT host_cpu QUERY PROCESSOR_DESCRIPTION)
math(EXPR host_ramGB "${host_ramMB} / 1000")
message(STATUS "${host_ramGB} GB RAM detected on ${CMAKE_HOST_SYSTEM_NAME} with ${host_cpu}.  Detected ${Ncpu} CPU cores.")
if(host_ramGB LESS 2)
  message(WARNING "Minimum RAM is about 2 GB--some tests or simulations may fail due to small memory (RAM)")
endif()


if(realbits EQUAL 32)
  message(VERBOSE " 32-bit real precision")
  set(arith s)
else()
  message(VERBOSE " 64-bit real precision")
  set(realbits 64)
  set(arith d)
endif()

option(dev "Gemini developer mode")

option(mpi "Use MPI parallelization" on)

option(autobuild "autobuild missing libraries" on)
option(glow "use NCAR GLOW airglow / aurora model" on)

option(hdf5_external "build HDF5 instead of finding")
option(lapack_external "build Lapack instead of finding")
option(scalapack_external "build ScaLapack instead of finding")
option(mumps_external "build MUMPS instead of finding")


option(hwm14 "use HWM14 neutral winds model")
option(msis20 "use MSIS 2.0 neutral atmosphere model")

option(hdf5 "use HDF5 file I/O" on)
option(netcdf "use NetCDF file I/O")
option(zlib_legacy "use unmaintained ZLIB 1.x")

# MUMPS build options (only used if auto-building MUMPS)
option(scotch "MUMPS: use Scotch" off)
option(openmp "MUMPS: use OpenMP" off)

option(python "PyGemini checks")
# Matlab checks take much longer than Python, and Python covers much more
option(matlab "Matlab checks")

set(CMAKE_TLS_VERIFY true)  # for Git and Downloads

# to make Gemini3D more usable by external programs, put all Fortran .mod generated module files in a single directory.
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

if(EXISTS ${PROJECT_SOURCE_DIR}/../mat_gemini/setup.m)
  set(FETCHCONTENT_SOURCE_DIR_MATGEMINI ${PROJECT_SOURCE_DIR}/../mat_gemini CACHE PATH "MatGemini developer path")
endif()

if(dev)

else()
  set_directory_properties(PROPERTIES EP_UPDATE_DISCONNECTED true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_MSIS2 true)
endif()

# --- auto-ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()

# --- default install directory under build/local
# users can specify like "cmake -B build -DCMAKE_INSTALL_PREFIX=~/mydir"
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # will not take effect without FORCE
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE PATH "Install top-level directory" FORCE)
endif()

# --- special handling of MacOS Homebrew/Macports
include(${CMAKE_CURRENT_LIST_DIR}/macos.cmake)
