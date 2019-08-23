mpimod = library('mpimod',
  sources: ['mpimod/mpimod.F90', 'mpimod/mpisend.f90', 'mpimod/mpirecv.f90', 'mpimod/mpihalo.f90'],
  dependencies: mpi,
  link_with: const,
  fortran_args: real_bits)

grid = library('grid',
  sources: ['grid/grid.f90', 'grid/read.f90'],
  dependencies: mpi,
  link_with: [const, mpimod])

calculus = library('calculus', 'calculus/calculus.f90',
  link_with: [const, grid])

advection = library('advection', 'advection/advec_mpi.f90',
  link_with: [const, grid, mpimod])

subdir('diffusion')
subdir('interpolation')
subdir('potential')