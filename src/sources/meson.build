sources = library('sources',
  sources: ['sources.f90', 'sources_mpi.f90'],
  link_with: [calculus, collisions, const, grid, mesh, mpimod],
    dependencies: mpi)  # ifort needs dep: mpi here