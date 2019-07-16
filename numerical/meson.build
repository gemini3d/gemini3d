mpi_mod = library('mpimod',
  'mpimod/mpimod.F90', 'mpimod/mpisend.f90', 'mpimod/mpirecv.f90', 'mpimod/mpihalo.f90',
  dependencies: mpi,
  fortran_args: REALBITS)

const = library('const', 'constants/phys_consts.F90')