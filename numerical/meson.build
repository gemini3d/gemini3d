mpimod = library('mpimod',
  'mpimod/mpimod.F90', 'mpimod/mpisend.f90', 'mpimod/mpirecv.f90', 'mpimod/mpihalo.f90',
  dependencies: mpi,
  link_with: const,
  fortran_args: real_bits)

grid = library('grid', 'grid/grid.f90',
  dependencies: mpi,
  link_with: [const, mpimod])

calculus = library('calculus', 'calculus/calculus.f90', 'interpolation/interpolation.f90',
  link_with: [const, grid, mpimod],
  dependencies : mpi,
  fortran_args: real_bits)

subdir('potential')