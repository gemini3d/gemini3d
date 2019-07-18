project('GEMINI', 'fortran',
  meson_version: '>=0.51.2',
  default_options : ['default_library=static', 'buildtype=release'])

REALBITS = '-DREALBITS='+get_option('realbits')
GLOW = get_option('GLOW')
HDF = get_option('HDF')
mpi_root = get_option('MPI_ROOT')
lapack_root = get_option('LAPACK_ROOT')
scalapack_root = get_option('SCALAPACK_ROOT')
mumps_root = get_option('MUMPS_ROOT')

fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_global_arguments('-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function', language : 'fortran')
  oldflags = ['-std=legacy', '-w']
else
  oldflags = ['-w']
endif

# need find_library to search specific dirs at this time
if mpi_root == ''
  mpi = dependency('mpi', language : 'fortran', required : true)
else
  message('searching for MPI in ' + mpi_root)
  mpi = fc.find_library('mpi', required : true, dirs : mpi_root / 'lib')
endif

mpiret = fc.run('use mpi; end', dependencies : mpi)
if mpiret.returncode() != 0
  error('MPI does not link properly: ' + mpiret.stderr() + mpiret.stdout() + mpi_root + ' ' + mpi.name() + ' ' + mpi.version())
endif

if lapack_root == ''
  lapack = dependency('lapack', required : true)
else
  message('searching for LAPACK in ' + lapack_root)
  lapack = fc.find_library('lapack', required : true, dirs : lapack_root)
endif

if HDF  # see HDF5 branch, not yet fully ready.
 hdf5 = dependency('hdf5')
endif

#==========================================================
const = library('const', 'numerical/constants/phys_consts.F90',
  fortran_args: REALBITS)

subdir('io')
subdir('temporal')
subdir('tests')

collisions = library('collisions', 'collisions/collisions.f90',
  link_with: const)

subdir('numerical')
io = library('io','io/io.F90',
  link_with: [calculus, const, fsutils, grid, mpimod])

msis = library('msis', 'vendor/msis00/msis00_gfortran.f',
  fortran_args : oldflags)

neutral = library('neutral', 'neutral/neutral.f90',
  dependencies : mpi,
  link_with : [const, calculus, mpimod, grid, io, msis, timeutils])

subdir('ionization')

# -- gemini exe
gemini = executable('gemini.bin', 'gemini.f90', 'temporal/temporal.f90', 'multifluid/multifluid.f90',
'sources/sources.f90', 'numerical/advection/advec_mpi.f90', 'numerical/diffusion/diffusion.F90', 'vendor/lapack95/gbsv.F90',
  link_with : [io,calculus, collisions, grid, ionization, mpimod, neutral, potential, const, timeutils],
  dependencies : [mpi, lapack],
  fortran_args: REALBITS)

# -- magcalc exe

magcalc = executable('magcalc.bin', 'magcalc.f90',
  link_with : [const, calculus, grid, neutral, io, timeutils, mpimod],
  dependencies : mpi)

# --- TESTS ---

octave = find_program('octave', required : false)
mpiexec = find_program('mpiexec', required : true)

test2d = meson.build_root() / 'test2d'
test3d = meson.build_root() / 'test3d'
zenodo2d_url = 'https://zenodo.org/record/2520780/files/zenodo2d.zip?download=1'
zenodo2d_md5 = '5bd1bce1a465ccec5af813f8b7959ec8'
zenodo2d_dir = meson.source_root() / '../simulations/zenodo2d'
zenodo2d_zip = zenodo2d_dir/'../zenodo2d.zip'
zenodo3d_url = 'https://zenodo.org/record/2520920/files/zenodo3d.zip?download=1'
zenodo3d_md5 = '225853d43937a70c9ef6726f90666645'
zenodo3d_dir = meson.source_root() / '../simulations/zenodo3d'
zenodo3d_zip = zenodo3d_dir/'../zenodo3d.zip'
eval2d = 'compare_all("' + test2d + '","' + zenodo2d_dir + '")'
eval3d = 'compare_all("' + test3d + '","' + zenodo3d_dir + '")'

run_command('python', '-c', 'import requests', check: true)

run_command('python', 'script_utils/meson_file_download.py',
  zenodo2d_url, zenodo2d_zip, zenodo2d_md5,
  check: true)
run_command('python', 'script_utils/meson_file_download.py',
  zenodo3d_url, zenodo3d_zip, zenodo3d_md5,
  check: true)

run_command('python', 'script_utils/meson_file_extract.py',
  zenodo2d_zip, zenodo2d_dir,
  check: true)
run_command('python', 'script_utils/meson_file_extract.py',
  zenodo3d_zip, zenodo3d_dir,
  check: true)

test('Gemini2D', mpiexec,
  suite : '2d',
  args : ['-np', '2', gemini, meson.source_root() / 'initialize/test2d/config.ini', test2d],
  workdir : meson.source_root(),
  timeout : 300,
  is_parallel : false)

if octave.found()
  test('Compare2D', octave, args : ['--eval', eval2d],
    suite : '2d', timeout: 30, workdir : meson.source_root() / 'tests')
endif

test('Gemini3D', mpiexec,
  suite : '3d',
  args : ['-np', '2', gemini, meson.source_root() / 'initialize/test3d/config.ini', test3d],
  workdir : meson.source_root(),
  timeout : 600,
  is_parallel : false)

if octave.found()
  test('Compare3D', octave, args : ['--eval', eval3d],
    suite : '3d', timeout: 30, workdir : meson.source_root() / 'tests')
endif
