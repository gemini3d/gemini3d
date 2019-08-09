project('GEMINI', 'fortran',
  meson_version: '>=0.51.2',
  default_options : ['default_library=static', 'buildtype=release',
  'libdir=gemini/lib', 'includedir=gemini/include'])

real_bits = '-DREALBITS=' + get_option('realbits')
use_glow = get_option('useglow')
use_hdf = get_option('HDF')

fc = meson.get_compiler('fortran')
# -- compiler options
if fc.get_id() == 'gcc'
  add_project_arguments('-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function',
  '-fimplicit-none',
  language : 'fortran')
  if fc.version().version_compare('>=8.0')
    add_project_arguments('-std=f2018', language: 'fortran')
  endif
elif fc.get_id() == 'intel'
  add_project_arguments('-stand', 'f18', '-warn declarations', language: 'fortran')
elif fc.get_id() == 'intel-cl'
   add_project_arguments('/stand:f18', '/warn:declarations', language: 'fortran')
endif
# need find_library to search specific dirs at this time
mpi_root = get_option('MPI_ROOT')
if mpi_root == ''
  mpi = dependency('mpi', language : 'fortran', required : true)
else
  message('searching for MPI in ' + mpi_root)
  mpi = fc.find_library('mpi', required : true, dirs : mpi_root / 'lib')
endif

mpiret = fc.run('use mpi; end', dependencies : mpi, name: 'Fortran MPI')
if mpiret.returncode() != 0
  error('MPI does not link properly: ' + mpiret.stderr() + mpiret.stdout() + mpi_root + ' ' + mpi.name() + ' ' + mpi.version())
endif

lapack_root = get_option('LAPACK_ROOT')
if lapack_root == ''
  lapack = dependency('lapack', required: false)
  blas = dependency('blas', required: false)
  if not blas.found()  # necessary for some systems incl. CentOS 7
    blas = fc.find_library('blas', required: false)
  endif
else
  message('searching for LAPACK & BLAS in ' + lapack_root)
  lapack = fc.find_library('lapack', dirs : lapack_root, required: false)
  blas = fc.find_library('blas', dirs : lapack_root, required: false)
endif
if not lapack.found() or not blas.found()
  lapack_proj = subproject('lapack')
  lapack = lapack_proj.get_variable('lapack')
  blas = lapack_proj.get_variable('blas')
endif

if use_hdf  # see HDF5 branch, not yet fully ready.
 hdf5 = dependency('hdf5')
endif

#==========================================================
const = library('const', 'numerical/constants/phys_consts.F90',
  fortran_args: real_bits)

subdir('io')
subdir('temporal')
subdir('tests')

collisions = library('collisions', 'collisions/collisions.f90',
  link_with: const)

subdir('numerical')
io = library('io','io/io.F90',
  link_with: [calculus, const, fsutils, grid, mpimod])

subdir('vendor/msis00')
# --- for setting up an equilibrium simulation --
executable('msis_setup', 'setup/MSIS00/call_msis_gfortran.f90',
  link_with: msis)

neutral = library('neutral', 'neutral/neutral.f90',
  dependencies : mpi,
  link_with : [const, calculus, mpimod, grid, io, msis, timeutils])

subdir('ionization')

gemini_src = files('gemini.f90', 'temporal/temporal.f90', 'multifluid/multifluid.f90',
'sources/sources.f90', 'numerical/advection/advec_mpi.f90', 'numerical/diffusion/diffusion.F90', 'vendor/lapack95/gbsv.F90')
# -- gemini exe
gemini_fang_exe = executable('gemini_fang.bin',
  sources: gemini_src,
  link_with : [io, potential, calculus, collisions, grid, ionization, mpimod, neutral, const, timeutils],
  dependencies : [scalapack, blacs, lapack, blas, mpi],
  fortran_args: real_bits,
  install: true,
  link_language: 'fortran')

# -- gemini_glow exe

if use_glow
  gemini_glow_exe = executable('gemini_glow.bin',
    sources: gemini_src,
    link_with : [io, potential, calculus, collisions, grid, ionization_glow, mpimod, neutral,  const, timeutils],
    dependencies : [scalapack, blacs, lapack, blas, mpi],
    fortran_args: real_bits,
    install: true,
    link_language: 'fortran')
else
  gemini_glow_exe = disabler()
endif

# -- magcalc exe

magcalc = executable('magcalc.bin', 'magcalc.f90',
  link_with : [const, calculus, grid, neutral, io, timeutils, mpimod],
  dependencies : mpi,
  install: true)

# --- TESTS ---
mpiexec = find_program('mpiexec')
pymod = import('python')
python = pymod.find_installation('python3')
octave = find_program('octave', required: false)

matlab = disabler()
if get_option('usematlab')  # buggy for Matlab < R2019a, may hang "meson setup" silently for two minutes
  matlab = find_program('matlab', required: false)
  if run_command(python, 'script_utils/meson_matlab_has_batch.py').returncode() != 0
    matlab = disabler()
  endif
endif

is_ci = run_command(python, 'script_utils/meson_get_environment_variable.py', 'CI').stdout().to_lower() == 'true'

zenodo = {
  '2d_fang': {
    'url': 'https://zenodo.org/record/2520780/files/zenodo2d.zip?download=1',
    'md5': '5bd1bce1a465ccec5af813f8b7959ec8',
    'dir': meson.source_root() / 'tests/data/zenodo2d',
    'zip': meson.source_root() / 'tests/data/zenodo2d.zip',
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

  if not use_glow and k.endswith('glow')
    continue
  endif

  if is_ci and k.startswith('3d')
    # save download / CPU on CI
    continue
  endif

  proc = run_command(python, files('script_utils/meson_file_download.py'),
                     v['url'], v['zip'], '-hash', 'md5', v['md5'])
  download_ok = true
  if proc.returncode() != 0
    warning(k + ' ' + v['zip'] + ' test data could not be downloaded from ' + v['url'] + ' and so ' + k + ' test will fail.')
    message('stderr: \n' + proc.stderr())
    download_ok = false
  endif

  ref_ok = true
  if download_ok
    proc = run_command(python, files('script_utils/meson_file_extract.py'),
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
    proc = run_command(python,
                       files('script_utils/meson_cpu_count.py'),
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
    args : ['-np', mpi_np, v['exe'], meson.source_root() / 'initialize/test' + k / 'config.ini', test_dir],
    workdir : meson.source_root(),
    timeout : 1200,
    is_parallel : false)

  if ref_ok
    test('Compare' + k, python,
      args : ['compare_all.py', test_dir, v['dir']],
      suite : k, timeout: 60, workdir : meson.source_root() / 'tests')

    matoctargs = 'compare_all("' + test_dir + '","' + v['dir'] + '")'

    if matlab.found()
      test('CompareMatlab' + k, matlab,
        args : ['-batch', matoctargs],
        suite : k, timeout: 60, workdir : meson.source_root() / 'tests')
    endif

    if octave.found()
      test('CompareOctave' + k, octave,
        args : ['--eval', matoctargs],
        suite : k, timeout: 60, workdir : meson.source_root() / 'tests')
    endif
  else
    warning('Compare' + k + ' skipped due to missing reference data in ' + v['dir'])
  endif

endforeach
