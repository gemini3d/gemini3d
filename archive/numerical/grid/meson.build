mesh = library('mesh',
  sources: 'mesh.f90',
  link_with: const)

grid_src = files('grid.f90', 'read.f90')
grid_deps = [mpi]
if hdf5_interface.found()
  grid_src += files('readgrid_hdf5.f90')
  grid_deps += hdf5_interface
elif netcdf_interface.found()
  grid_src += files('readgrid_nc4.f90')
  grid_deps += netcdf_interface
else
  grid_src += files('readgrid_raw.f90')
endif

grid = library('grid',
  sources: grid_src,
  dependencies: grid_deps,
  link_with: [const, reader, mesh, mpimod])