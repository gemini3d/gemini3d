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

cmake_host_system_information(RESULT _ramMiB QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT _cpu QUERY PROCESSOR_DESCRIPTION)
math(EXPR _ramGB "${_ramMiB} / 1000")
message(STATUS "${_ramGB} GB RAM detected on ${CMAKE_HOST_SYSTEM_NAME} with ${_cpu}.  Detected ${Ncpu} CPU cores.")
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

# If MPIEXEC is not present, no need to link MPI.
if(NOT DEFINED mpi)
  find_program(_mpiexec NAMES mpiexec)
  if(_mpiexec)
    set(mpi on)
  else()
    set(mpi off)
  endif(_mpiexec)
endif()

option(dev "Gemini developer mode")

option(mpi "Use MPI parallelization")

option(autobuild "autobuild missing Lapack, Scalapack or Mumps" on)

option(hdf5_external "build HDF5 instead of finding")
option(lapack_external "build Lapack instead of finding")
option(scalapack_external "build ScaLapack instead of finding")
option(mumps_external "build MUMPS instead of finding")

option(glow "use NCAR GLOW airglow / aurora model" on)
option(hwm14 "use HWM14 neutral winds model")
option(msis20 "use MSIS 2.0 neutral atmosphere model")

option(hdf5 "use HDF5 file I/O" on)
option(netcdf "use NetCDF file I/O")

# MUMPS build options (only used if auto-building MUMPS)
option(metis "MUMPS: use METIS" off)
option(scotch "MUMPS: use Scotch" off)
option(openmp "MUMPS: use OpenMP" off)

option(python "PyGemini checks")
# Matlab checks take much longer than Python, and Python covers much more
option(matlab "Matlab checks")

if(dev)
  set(FETCHCONTENT_SOURCE_DIR_PYGEMINI ${PROJECT_SOURCE_DIR}/../pygemini CACHE PATH "PyGemini developer path")
  set(FETCHCONTENT_SOURCE_DIR_MATGEMINI ${PROJECT_SOURCE_DIR}/../mat_gemini CACHE PATH "MatGemini developer path")
else()
  set(FETCHCONTENT_UPDATES_DISCONNECTED_GLOW true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_H5FORTRAN true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_HWM14 true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_LAPACK true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_MATGEMINI true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_MSIS2 true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_MUMPS true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_NC4FORTRAN true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_PYGEMINI true)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_SCALAPACK true)
endif()
