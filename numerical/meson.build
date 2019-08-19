mpimod = library('mpimod',
  sources: ['mpimod/mpimod.F90', 'mpimod/mpisend.f90', 'mpimod/mpirecv.f90', 'mpimod/mpihalo.f90'],
  dependencies: mpi,
  link_with: const,
  fortran_args: real_bits)

grid = library('grid', 'grid/grid.f90',
  dependencies: mpi,
  link_with: [const, mpimod])

calculus = library('calculus', 'calculus/calculus.f90',
  link_with: [const, grid])

advection = library('advection', 'advection/advec_mpi.f90',
  link_with: [const, grid, mpimod])

diffusion = library('diffusion',
  sources: ['diffusion/diffusion.F90', '../vendor/lapack95/gbsv.F90'],
  dependencies: lapack,
  link_with: [const, grid],
  fortran_args: real_bits)

interp = library('interpolation', 'interpolation/interpolation.f90',
  link_with: const)

subdir('potential')