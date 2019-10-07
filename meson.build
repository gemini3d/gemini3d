project('GEMINI', 'fortran',
  meson_version: '>=0.52.0',
  default_options : ['default_library=static', 'buildtype=release',
  'libdir=gemini/lib', 'includedir=gemini/include'])

bits = get_option('realbits')
real_bits = '-DREALBITS=' + bits
use_hdf = get_option('HDF')
buildtype = get_option('buildtype')
arith = get_option('arith')
message('arith: ' + arith)

subdir('cmake')

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
  link_with: [const, grid, mpimod])

subdir('tests')

io = library('io',
  sources: io_src,
  link_with: [calculus, const, fsutils, grid, mpimod, timeutils])

subdir('vendor/msis00')
# --- for setting up an equilibrium simulation --
executable('msis_setup', 'setup/MSIS00/call_msis_gfortran.f90',
  link_with: msis)

neutral = library('neutral', 'neutral/neutral.f90',
  dependencies : mpi,
  link_with : [const, grid, interp, mpimod, msis, timeutils])

subdir('ionization')

sources = library('sources', 'sources/sources.f90',
  link_with: [calculus, collisions, const, grid, mpimod])

multifluid = library('multifluid', 'multifluid/multifluid.f90',
  link_with: [advection, calculus, collisions, const, diffusion, grid, ionization_fang, mpimod,
              precipBCs, sources, timeutils])


# -- gemini exe
gemini_fang_exe = executable('gemini_fang.bin',
  sources: 'gemini.f90',
  link_with : [const, grid, io, mpimod, multifluid, neutral, potential, precipBCs, temporal, timeutils],
  dependencies : [scalapack, blacs, lapack, blas],
  fortran_args: real_bits,
  install: true,
  link_language: 'fortran')

# -- gemini_glow exe

multifluid_glow = library('multifluid_glow', 'multifluid/multifluid.f90',
  link_with: [advection, calculus, collisions, const, diffusion, grid, ionization_glow, mpimod,
              precipBCs, sources, timeutils])

gemini_glow_exe = executable('gemini_glow.bin',
  sources: 'gemini.f90',
  link_with : [const, grid, io, mpimod, multifluid_glow, neutral, potential, precipBCs, temporal, timeutils],
  dependencies : [scalapack, blacs, lapack, blas],
  fortran_args: real_bits,
  install: true,
  link_language: 'fortran')

# -- magcalc exe

magcalc = executable('magcalc.bin', 'magcalc.f90',
  link_with : [const, calculus, grid, neutral, io, timeutils, mpimod],
  dependencies : mpi,
  install: true)

# --- TESTS ---

zenodo = {
  '2d_fang': {
    'url': 'https://zenodo.org/record/3464571/files/zenodo2d_fang.zip?download=1',
    'md5': '57d72fd0005247c8eedf122ac4670ad0',
    'dir': meson.source_root() / 'tests/data/zenodo2d_fang',
    'zip': meson.source_root() / 'tests/data/zenodo2d_fang.zip',
    'exe': gemini_fang_exe},
  '2d_glow': {
    'url': 'https://zenodo.org/record/2520780/files/zenodo2d_glow.zip?download=1',
    'md5': 'c5bbbbff3bdde85b6d7e9470bc3751a2',
    'dir': meson.source_root() / 'tests/data/zenodo2d_glow',
    'zip': meson.source_root() / 'tests/data/zenodo2d_glow.zip',
    'exe': gemini_glow_exe},
  '3d_fang': {
    'url': 'https://zenodo.org/record/2520920/files/zenodo3d.zip?download=1',
    'md5': '225853d43937a70c9ef6726f90666645',
    'dir': meson.source_root() / 'tests/data/zenodo3d',
    'zip': meson.source_root() / 'tests/data/zenodo3d.zip',
    'exe': gemini_fang_exe},
  '3d_glow': {
     'url': 'https://zenodo.org/record/2520920/files/zenodo3d_glow.zip?download=1',
     'md5': '2c54bfde8aff0fb72d61115f04c361a7',
     'dir': meson.source_root() / 'tests/data/zenodo3d_glow',
     'zip': meson.source_root() / 'tests/data/zenodo3d_glow.zip',
     'exe': gemini_glow_exe},
  }

mpi_np = get_option('NP')

foreach k, v : zenodo

#  if is_ci and k.startswith('3d')
#    # save download / CPU on CI
#    continue
#  endif

  proc = run_command(python, meson.source_root() / 'script_utils/meson_file_download.py',
                     v['url'], v['zip'], '-hash', 'md5', v['md5'])
  download_ok = true
  if proc.returncode() != 0
    warning(k + ' ' + v['zip'] + ' test data could not be downloaded from ' + v['url'] + ' and so ' + k + ' test will fail.')
    message('stderr: \n' + proc.stderr())
    download_ok = false
  endif

  ref_ok = true
  if download_ok
    proc = run_command(python, meson.source_root() / 'script_utils/meson_file_extract.py',
                       v['zip'], v['dir'])
    if proc.returncode() != 0
      warning(k + 'test data in ' + v['dir'] + ' not available: ' + k + ' test disabled.')
      message('stderr: \n' + proc.stderr())
      message('stdout: \n' + proc.stdout())
      ref_ok=false
    endif
  endif

  test_dir = meson.build_root() / 'test'  + k

  if mpi_np == ''
    proc = run_command(python, meson.source_root() / 'script_utils/meson_cpu_count.py',
                       v['dir']/ 'inputs/simsize.dat')
    mpi_np = proc.stdout()
    if proc.returncode() != 0 or mpi_np == ''
      warning('Unable to auto-compute appropriate number of MPI images, using NP=2')
      mpi_np = '2'
    endif
  endif

  testname = 'Gemini' + k + '-NP' + mpi_np
  message('Enabled test ' + testname)
  test(testname, mpiexec,
    suite : k,
    args : ['-np', mpi_np, v['exe'], meson.source_root() / 'initialize/test' + k / 'config.nml', test_dir],
    workdir : meson.source_root(),
    timeout : 1200,
    is_parallel : false)

  if ref_ok
    test('Compare' + k, python_test,
      args : ['tests/compare_all.py', test_dir, v['dir']],
      suite : k,
      timeout: 60,
      workdir : meson.source_root())

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
  else
    warning('Compare' + k + ' skipped due to missing reference data in ' + v['dir'])
  endif

endforeach
