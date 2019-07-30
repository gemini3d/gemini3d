ionization_src = files('ionization.f90', 'boundary_conditions/precipBCs_mod.f90')

ionization = library('ionization', ionization_src, 'glow_dummy.f90',
  link_with: [io, const, calculus, grid, neutral, timeutils, mpimod],
  dependencies : mpi)

# ===== optional GLOW library

if use_glow
  glow_proj = subproject('cglow')
  glow_lib = glow_proj.get_variable('cglow')

  data_dir = 'subprojects/glow/data/' # need trailing slash

  ionization_glow = library('ionization_glow', ionization_src, 'glow_run.F90',
    link_with: [glow_lib, io, const, calculus, grid, neutral, timeutils, mpimod],
    dependencies : mpi,
    fortran_args: '-DDATADIR="'+data_dir+'"')
endif