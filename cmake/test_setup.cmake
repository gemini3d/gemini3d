include(${CMAKE_CURRENT_LIST_DIR}/test_compare.cmake)


function(get_mpi_count outdir)
# not very optimal, use test runner instead (Fortran, Python or Matlab)
find_program(h5dump NAMES h5dump)
if(NOT h5dump)
  message(STATUS "skipping test ${testname} because h5dump is not found")
  return()
endif()

set(_tmpfile ${PROJECT_BINARY_DIR}/lx.txt)

execute_process(COMMAND ${h5dump} -o ${_tmpfile} --noindex ${outdir}/inputs/simsize.h5
  OUTPUT_QUIET
  RESULT_VARIABLE _err
  TIMEOUT 10)
if(NOT _err EQUAL 0)
  message(STATUS "SKIP: ${_err} ${testname} due to problem reading ${outdir}/inputs/simsize.h5")
  return()
endif()

file(STRINGS ${_tmpfile} lx REGEX "[0-9]+" LIMIT_COUNT 3 LIMIT_INPUT 80)
list(GET lx 1 lx2)
list(GET lx 2 lx3)

if(lx2 GREATER lx3)
  set(_lxm ${lx2})
else()
  set(_lxm ${lx3})
endif()

if(${MPIEXEC_MAX_NUMPROCS} GREATER ${_lxm})
  set(mpi_count ${_lxm})
else()
  math(EXPR _c "${_lxm} % ${MPIEXEC_MAX_NUMPROCS}")
  if(_c EQUAL 0)
    set(mpi_count ${_lxm})
  else()
    math(EXPR mpi_count "${MPIEXEC_MAX_NUMPROCS} - ${_c}")
  endif()

  if(mpi_count LESS 1)
    set(mpi_count 1)
  endif()
endif()

message(STATUS "${testname} MPI images: ${mpi_count}")

set(mpi_count ${mpi_count} PARENT_SCOPE)

endfunction(get_mpi_count)


function(setup_gemini_test testname TIMEOUT)

# --- setup test
set(_outdir ${PROJECT_BINARY_DIR}/test${testname})
set(_refdir ${PROJECT_SOURCE_DIR}/tests/data)

if(PYGEMINI_DIR)

  if(mpi)
    set(_cmd ${Python3_EXECUTABLE} -m gemini3d.run_test -mpiexec ${MPIEXEC_EXECUTABLE} ${testname} $<TARGET_FILE:gemini.bin> ${_outdir} ${_refdir})
  else()
    set(_cmd ${Python3_EXECUTABLE} -m gemini3d.run_test ${testname} $<TARGET_FILE:gemini.bin> ${_outdir} ${_refdir})
  endif(mpi)

else(PYGEMINI_DIR)

  set(_r ${_refdir}/test${testname}/inputs)

  if(NOT IS_DIRECTORY ${_outdir})
    if(IS_DIRECTORY ${_r})
      file(COPY ${_r} DESTINATION ${PROJECT_BINARY_DIR}/test${testname})
    else()
      # message(STATUS "SKIP: ${testname} data not in ${_refdir}. Run cmake -Ddownload=yes to get ref data.")
      # return()
    endif()
  endif()

  if(mpi)
    set(_cmd $<TARGET_FILE:gemini3d.run> ${_outdir} -n ${MPIEXEC_MAX_NUMPROCS} -gemexe $<TARGET_FILE:gemini.bin> -mpiexe ${MPIEXEC_EXECUTABLE})
  else()
    set(_cmd $<TARGET_FILE:gemini.bin> ${_outdir})
  endif(mpi)

endif(PYGEMINI_DIR)

if(hdf5)

add_test(NAME gemini:hdf5:${testname}:dryrun
  COMMAND ${_cmd} -dryrun)
# we prefer default WorkingDirectory of PROJECT_BINARY_DIR to make MSIS 2.0 msis20.parm use simpler
# otherwise, we have to generate source for msis_interface.f90

  # WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  # NOTE: WorkingDirectory is NECESSARY for Windows + Intel + ProgramFiles/HDFGroup/HDF5

set_tests_properties(gemini:hdf5:${testname}:dryrun PROPERTIES
  TIMEOUT 60
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED mumps_fixture
  FIXTURES_SETUP hdf5:${testname}:dryrun)


add_test(NAME gemini:hdf5:${testname}
  COMMAND ${_cmd})
  # WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

# NOTE: don't use REQUIRED_FILES because it won't let file download if not present.
set_tests_properties(gemini:hdf5:${testname} PROPERTIES
  TIMEOUT ${TIMEOUT}
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED hdf5:${testname}:dryrun
  FIXTURES_SETUP hdf5:${testname})

endif(hdf5)

if(netcdf)
add_test(NAME gemini:netcdf:${testname}:dryrun
  COMMAND ${_cmd} -out_format nc -dryrun)

set_tests_properties(gemini:netcdf:${testname}:dryrun PROPERTIES
  TIMEOUT 60
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED mumps_fixture
  FIXTURES_SETUP netcdf:${testname}:dryrun)

add_test(NAME gemini:netcdf:${testname}
  COMMAND ${_cmd} -out_format nc)

set_tests_properties(gemini:netcdf:${testname} PROPERTIES
  TIMEOUT ${TIMEOUT}
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED netcdf:${testname}:dryrun
  FIXTURES_SETUP netcdf:${testname})
endif(netcdf)

compare_gemini_output(${testname})

endfunction(setup_gemini_test)
