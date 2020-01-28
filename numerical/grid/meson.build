mesh = library('mesh', 'mesh.f90',
  link_with: const)

grid = library('grid',
  sources: ['grid.f90', 'read.f90'],
  dependencies: mpi,
  link_with: [const, reader, mesh, mpimod])