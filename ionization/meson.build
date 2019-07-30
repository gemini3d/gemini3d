ionization_src = files('ionization.f90', 'boundary_conditions/precipBCs_mod.f90')

if use_glow
  ionization_src += 'glow_run.F90'
  glow_proj = subproject('cglow')
  cglow = glow_proj.get_variable('cglow')
  ionization_deps = cglow
else
  ionization_src += 'glow_dummy.f90'
  ionization_deps = []
endif

data_dir = 'subprojects/glow/data/' # need trailing slash

ionization = library('ionization', ionization_src,
  link_with: [ionization_deps, io, const, calculus, grid, neutral, timeutils, mpimod],
  dependencies : mpi,
  fortran_args: '-DDATADIR="'+data_dir+'"')