fang = library('ionize_fang', 'fang.f90',
  link_with: const)

test_fang2008_exe = executable('test_fang2008', 'test_fang2008.f90',
  link_with: [fang, msis])

test('Fang2008', python,
  args: ['test_fang2008.py', test_fang2008_exe],
  timeout: 10,
  workdir: meson.current_source_dir(),
  suite: 'unit')

ionization_fang = library('ionization',
  sources: ['ionization.f90', 'glow_dummy.f90'],
  link_with: [const, neutral, grid, timeutils, mpimod, fang],
  dependencies : mpi)

precipBCs = library('precipBCs', 'boundary_conditions/precipBCs_mod.f90',
  link_with: [mpimod, const, grid, interp, timeutils],
  dependencies: mpi)


# ===== GLOW library

glow_proj = subproject('cglow')
glow_lib = glow_proj.get_variable('cglow')

data_dir = 'subprojects/glow/data/' # need trailing slash

ionization_glow = library('ionization_glow',
  sources: ['ionization.f90', 'glow_run.F90'],
  link_with: [glow_lib, const, neutral, grid, timeutils, mpimod, fang],
  dependencies : mpi,
  fortran_args: '-DDATADIR="'+data_dir+'"')
