project('Gemini3D', 'fortran',
  version: '0.3.1',
  meson_version: '>=0.52.0',
  default_options : ['default_library=static', 'buildtype=release'])

message('CMake is recommended instead of Meson, which is not maintained by PyGemini.')

subdir('cmake/meson')  # find libraries, etc.

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

io_deps = [mpi]
if hdf5_interface.found()
  io_deps += hdf5_interface
elif netcdf_interface.found()
  io_deps += netcdf_interface
endif
io = library('io',
  sources: io_src,
  link_with: [calculus, const, pathlib, reader, grid, mpimod, timeutils],
  dependencies: io_deps)  # ifort needs dep: mpi here


subdir('vendor/msis00')
# --- for setting up an equilibrium simulation --
executable('msis_setup', 'setup/MSIS00/msis_driver.f90',
  link_with: msis)

subdir('neutral')
subdir('ionization')
subdir('sources')

multifluid = library('multifluid', 'multifluid/multifluid.f90',
  link_with: [advection, calculus, collisions, const, diffusion, grid, mesh, ionization, mpimod,
              precipBCs, sources, timeutils],
  dependencies: mpi)  # ifort needs dep: mpi here


# -- gemini exe
gemini_exe = executable('gemini.bin',
  sources: 'gemini.f90',
  link_with : [const, grid, mesh, io, mpimod, multifluid, neutral, potential, precipBCs, temporal, timeutils],
  dependencies : [scalapack, blacs, lapack, blas, mpi],  # ifort needs dep: mpi here
  install: true,
  link_language: 'fortran')

# -- magcalc exe

magcalc = executable('magcalc.bin', 'magcalc.f90',
  link_with : [const, calculus, grid, mesh, neutral, io, timeutils, mpimod],
  dependencies : mpi,
  install: true)

# --- TESTS ---

zenodo = {
  '2d_fang': meson.source_root() / 'tests/data/test2d_fang',
  '2d_glow': meson.source_root() / 'tests/data/test2d_glow',
  '3d_fang': meson.source_root() / 'tests/data/test3d_fang',
  '3d_glow': meson.source_root() / 'tests/data/test3d_glow'}

np = get_option('np')

foreach k, v : zenodo

  test_dir = meson.build_root() / 'test'  + k

  testname = 'Gemini' + k
  test_args = ['meson_run_test.py', k, mpiexec.path(), gemini_exe,
               meson.source_root() / 'initialize/test' + k / 'config.nml',
               test_dir]
  if np > 0
    test_args += ['-np', np.to_string()]
  endif

  test(testname, python,
    suite : k,
    args : test_args,
    workdir : meson.source_root() / 'script_utils',
    timeout : 2700,
    env: nomalloc_env,
    is_parallel : false)


  # --- Python-based simulation output comparision
  test('Compare' + k, python,
    args : ['compare_all.py', test_dir, v],
    suite : k,
    timeout: 60,
    env: nomalloc_env,
    workdir : meson.source_root() / 'script_utils')
  # --- Matlab/Octave-based simulation output comparision
  matoctargs = 'compare_all("' + test_dir + '","' + v + '")'

  test('CompareMatlab' + k, matlab,
    args : ['-batch', matoctargs],
    suite : k,
    timeout: 60,
    env: nomalloc_env,
    workdir : meson.source_root() / 'script_utils')

  test('CompareOctave' + k, octave,
    args : ['--eval', matoctargs],
    suite : k,
    timeout: 60,
    env: nomalloc_env,
    workdir : meson.source_root() / 'script_utils')

endforeach
