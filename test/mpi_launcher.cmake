function(test_mpi_command Nworker working_dir out_var)
# Can't use TEST_LAUNCHER or CROSSCOMPILING_EMULATOR because multiple tests with different Nworker overwrite the property for other tests.

if(NOT MPIEXEC_EXECUTABLE OR NOT MPIEXEC_NUMPROC_FLAG)
  message(FATAL_ERROR "MPIEXEC_EXECUTABLE and MPIEXEC_NUMPROC_FLAG are required to run MPI tests.")
endif()

if(NOT Nworker GREATER 0)
  message(FATAL_ERROR "Number of MPI workers must be strictly positive integer")
endif()

set(_mpi_cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Nworker})

# cannot be IS_DIRECTORY because working_dir may be a generator expression
if(NOT "${working_dir}" STREQUAL "")
  list(APPEND _mpi_cmd -wdir ${working_dir})
  message(DEBUG "Working directory ${working_dir} for test ${test}.")
endif()

set(${out_var} ${_mpi_cmd} PARENT_SCOPE)

endfunction()


function(test_mpi_props test Nworker)

set_property(TEST ${test} PROPERTY PROCESSORS ${Nworker})

if(DEFINED mpi_tmpdir)
  set_property(TEST ${test} PROPERTY ENVIRONMENT "TMPDIR=${mpi_tmpdir}")
endif()

endfunction()


function(test_mpi_launcher target test Nworker)

if(ARGC GREATER 3)
  set(working_dir ${ARGV3})
else()
  set(working_dir "")
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
  get_property(_launcher TARGET ${target} PROPERTY TEST_LAUNCHER)
else()
  get_property(_launcher TARGET ${target} PROPERTY CROSSCOMPILING_EMULATOR)
endif()
if(_launcher)
  message(FATAL_ERROR "MPI launcher is already set for target ${target}. Cannot set it again.")
endif()

test_mpi_command(${Nworker} "${working_dir}" mpi_cmd)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
  set_property(TARGET ${target} PROPERTY TEST_LAUNCHER ${mpi_cmd})
else()
  set_property(TARGET ${target} PROPERTY CROSSCOMPILING_EMULATOR ${mpi_cmd})
endif()

test_mpi_props(${test} ${Nworker})

endfunction()
