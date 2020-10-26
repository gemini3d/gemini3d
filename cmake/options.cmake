cmake_host_system_information(RESULT _ramMiB QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT _cpu QUERY PROCESSOR_DESCRIPTION)
math(EXPR _ramGB "${_ramMiB} / 1024")
message(STATUS "${_ramGB} GB RAM detected on ${CMAKE_HOST_SYSTEM_NAME} with ${_cpu}")
if(_ramGB LESS 2)
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

option(mpi "Use MPI parallelization" on)

option(autobuild "autobuild missing Lapack, Scalapack or Mumps" on)
option(glow "use NCAR GLOW airglow / aurora model" on)
option(hdf5 "use HDF5 file I/O" on)
option(netcdf "use NetCDF file I/O" off)
# MUMPS build options (only used if auto-building MUMPS)
option(metis "MUMPS: use METIS" off)
option(scotch "MUMPS: use Scotch" off)
option(openmp "MUMPS: use OpenMP" off)

# on: debug, off: normal
set(FETCHCONTENT_UPDATES_DISCONNECTED off)

# this helps linters e.g. Visual Studio Intellicode work properly
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
