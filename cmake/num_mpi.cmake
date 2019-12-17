function(num_mpi_processes REFDIR)

find_package(Python3 COMPONENTS Interpreter)

# do not quote COMMAND line in general
execute_process(
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/script_utils/meson_cpu_count.py ${REFDIR}/inputs/simsize.dat
  OUTPUT_VARIABLE NP
  TIMEOUT 15)

set(NP ${NP} PARENT_SCOPE)

endfunction(num_mpi_processes)

# for testing standalone cmake -P
# num_mpi_processes(tests/data/zenodo2d_fang)