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
message(STATUS "Gemini3D: ${host_ramGB} GB RAM detected on ${CMAKE_HOST_SYSTEM_NAME} with ${host_cpu}.  Detected ${Ncpu} CPU cores.")
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

option(cpp "also build Gemini3D C++ frontend prototype" on)

option(mpi "Use MPI parallelization" on)

option(glow "use NCAR GLOW airglow / aurora model" on)

option(hwm14 "use HWM14 neutral winds model")
option(msis2 "use MSIS 2.x neutral atmosphere model")

option(hdf5 "use HDF5 file I/O" on)
option(netcdf "use NetCDF file I/O")

option(python "PyGemini checks")
# Matlab checks take much longer than Python, and Python covers much more
option(matlab "Matlab checks")

set(CMAKE_TLS_VERIFY true)  # for Git and Downloads

cmake_path(SET CMAKE_MODULE_PATH NORMALIZE ${CMAKE_CURRENT_LIST_DIR}/Modules)

# append .debug to debug libraries, because the computation speed penalty is so great
set(CMAKE_DEBUG_POSTFIX .debug)

# to make Gemini3D more usable by external programs, put all Fortran .mod generated module files in a single directory.
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# --- External project generator
if(CMAKE_GENERATOR STREQUAL "Ninja Multi-Config")
  set(EXTPROJ_GENERATOR "Ninja")
else()
  set(EXTPROJ_GENERATOR ${CMAKE_GENERATOR})
endif()

# --- auto-ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()

# --- default install directory under build/local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # will not take effect without FORCE
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE PATH "Install top-level directory" FORCE)
endif()
