mpimod = library('mpimod',
  sources: ['mpimod/mpimod.F90', 'mpimod/mpisend.f90', 'mpimod/mpirecv.f90', 'mpimod/mpihalo.f90'],
  dependencies: mpi,
  link_with: const,
  fortran_args: real_bits)

subdir('grid')

calculus = library('calculus',
  sources: ['calculus/calculus.f90', 'calculus/gradient.f90', 'calculus/integral.f90', 'calculus/div.f90'],
  link_with: [const, mesh])

advection = library('advection', 'advection/advec_mpi.f90',
  link_with: [const, grid, mesh, mpimod],
  dependencies: mpi)  # ifort needs dep: mpi here

subdir('diffusion')
subdir('interpolation')
subdir('potential')
