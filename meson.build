project('GEMINI', 'fortran',
  meson_version: '>=0.52.0',
  default_options : ['default_library=static', 'buildtype=release'])

subdir('meson')  # find libraries, etc.

#==========================================================
const = library('const', 'numerical/constants/phys_consts.F90',
  fortran_args: real_bits)

subdir('io')

collisions = library('collisions', 'collisions/collisions.f90',
  link_with: const)
subdir('temporal')
subdir('numerical')

temporal = library('temporal', 'temporal/temporal.f90',
  dependencies: mpi,
  link_with: [const, mesh, mpimod])

subdir('tests')

if hdf5_interface.found()
  io = library('io',
    sources: io_src,
    link_with: [calculus, const, fsutils, grid, mpimod, timeutils],
    dependencies: [mpi, hdf5_interface])  # ifort needs dep: mpi here
else
  io = library('io',
    sources: io_src,
    link_with: [calculus, const, fsutils, grid, mpimod, timeutils],
    dependencies: mpi)  # ifort needs dep: mpi here
endif

subdir('vendor/msis00')
# --- for setting up an equilibrium simulation --
executable('msis_setup', 'setup/MSIS00/call_msis_gfortran.f90',
  link_with: msis)

neutral = library('neutral', 'neutral/neutral.f90',
  dependencies : mpi,
  link_with : [const, grid, interp, mpimod, msis, timeutils])

subdir('ionization')

sources = library('sources', 'sources/sources.f90',
  link_with: [calculus, collisions, const, grid, mpimod],
    dependencies: mpi)  # ifort needs dep: mpi here

multifluid = library('multifluid', 'multifluid/multifluid.f90',
  link_with: [advection, calculus, collisions, const, diffusion, grid, ionization_fang, mpimod,
              precipBCs, sources, timeutils],
  dependencies: mpi)  # ifort needs dep: mpi here


# -- gemini exe
gemini_fang_exe = executable('gemini_fang.bin',
  sources: 'gemini.f90',
  link_with : [const, grid, io, mpimod, multifluid, neutral, potential, precipBCs, temporal, timeutils],
  dependencies : [scalapack, blacs, lapack, blas, mpi],  # ifort needs dep: mpi here
  fortran_args: real_bits,
  install: true,
  link_language: 'fortran')

# -- gemini_glow exe

multifluid_glow = library('multifluid_glow', 'multifluid/multifluid.f90',
  link_with: [advection, calculus, collisions, const, diffusion, grid, ionization_glow, mpimod,
              precipBCs, sources, timeutils],
  dependencies: mpi)  # ifort needs dep: mpi here

gemini_glow_exe = executable('gemini_glow.bin',
  sources: 'gemini.f90',
  link_with : [const, grid, io, mpimod, multifluid_glow, neutral, potential, precipBCs, temporal, timeutils],
  dependencies : [scalapack, blacs, lapack, blas, mpi],    # ifort needs dep: mpi here
  fortran_args: real_bits,
  install: true,
  link_language: 'fortran')

# -- magcalc exe

magcalc = executable('magcalc.bin', 'magcalc.f90',
  link_with : [const, calculus, grid, neutral, io, timeutils, mpimod],
  dependencies : mpi,
  install: true)

# --- TESTS ---

# tests can't be run without mpiexec.
if not mpiexec.found()
  subdir_done()
endif

zenodo = {
  '2d_fang': {
    'dir': meson.source_root() / 'tests/data/zenodo2d_fang',
    'exe': gemini_fang_exe},
  '2d_glow': {
    'dir': meson.source_root() / 'tests/data/zenodo2d_glow',
    'exe': gemini_glow_exe},
  '3d_fang': {
    'dir': meson.source_root() / 'tests/data/zenodo3d_fang',
    'exe': gemini_fang_exe},
  '3d_glow': {
     'dir': meson.source_root() / 'tests/data/zenodo3d_glow',
     'exe': gemini_glow_exe},
  }

foreach k, v : zenodo

  test_dir = meson.build_root() / 'test'  + k

  testname = 'Gemini' + k
  test(testname, python,
    suite : k,
    args : ['meson_run_test.py', k, mpiexec.path(), v['exe'], meson.source_root() / 'initialize/test' + k / 'config.nml', test_dir],
    workdir : meson.source_root() / 'script_utils',
    timeout : 2700,
    is_parallel : false)


  # --- Python-based simulation output comparision
  test('Compare' + k, python,
    args : ['compare_all.py', test_dir, v['dir']],
    suite : k,
    timeout: 60,
    workdir : meson.source_root() / 'tests')
  # --- Matlab/Octave-based simulation output comparision
  matoctargs = 'compare_all("' + test_dir + '","' + v['dir'] + '")'

  test('CompareMatlab' + k, matlab,
    args : ['-batch', matoctargs],
    suite : k,
    timeout: 60,
    workdir : meson.source_root() / 'tests')

  test('CompareOctave' + k, octave,
    args : ['--eval', matoctargs],
    suite : k,
    timeout: 60,
    workdir : meson.source_root() / 'tests')

endforeach
