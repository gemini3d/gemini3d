project('GEMINI', 'fortran',
  meson_version: '>=0.51.2',
  default_options : ['default_library=static', 'buildtype=release'])

REALBITS = '-DREALBITS='+get_option('realbits')
use_glow = get_option('useglow')
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
matlab = find_program('matlab', required : false)
if matlab.found()
  octmat = matlab
  evalbatch = '-batch'
elif octave.found()
  octmat = octave
  evalbatch = '--eval'
else
  octmat = []
  evalbatch = ''
  message('Neither Matlab or Octave were found, "Compare" tests will fail.')
endif
mpiexec = find_program('mpiexec', required : true)

zenodo = {
  '2d': {'url': 'https://zenodo.org/record/2520780/files/zenodo2d.zip?download=1',
         'md5': '5bd1bce1a465ccec5af813f8b7959ec8',
         'dir': meson.source_root() / '../simulations/zenodo2d',
         'zip': meson.source_root() / '../simulations/zenodo2d.zip'},
  '2d_glow': {'url': 'https://zenodo.org/record/2520780/files/zenodo2d_glow.zip?download=1',
              'md5': 'c5bbbbff3bdde85b6d7e9470bc3751a2',
              'dir': meson.source_root() / '../simulations/zenodo2d_glow',
              'zip': meson.source_root() / '../simulations/zenodo2d_glow.zip'},
  '3d': {'url': 'https://zenodo.org/record/2520920/files/zenodo3d.zip?download=1',
         'md5': '225853d43937a70c9ef6726f90666645',
         'dir': meson.source_root() / '../simulations/zenodo3d',
         'zip': meson.source_root() / '../simulations/zenodo3d.zip'},
  '3d_glow': {'url': 'https://zenodo.org/record/2520920/files/zenodo3d_glow.zip?download=1',
              'md5': '2c54bfde8aff0fb72d61115f04c361a7',
              'dir': meson.source_root() / '../simulations/zenodo3d_glow',
              'zip': meson.source_root() / '../simulations/zenodo3d_glow.zip'},
  }

if use_glow
  test_names = ['2d_glow', '3d_glow']
else
  test_names = ['2d', '3d']
endif


foreach i : test_names

  run_command('python', 'script_utils/meson_file_download.py',
    zenodo[i]['url'], zenodo[i]['zip'], 'md5', zenodo[i]['md5'])

  if run_command('python', 'script_utils/meson_file_extract.py',
    zenodo[i]['zip'], zenodo[i]['dir']).returncode() != 0
    message(i + 'test data not available, test will fail')
  endif

  test_dir = meson.build_root() / 'test'  + i

  test('Gemini' + i, mpiexec,
    suite : i,
    args : ['-np', '2', gemini, meson.source_root() / 'initialize/test' + i / 'config.ini', test_dir],
    workdir : meson.source_root(),
    timeout : 300,
    is_parallel : false)

  if octmat.found()
    octmat_args = [evalbatch, 'compare_all("' + test_dir + '","' + zenodo[i]['dir'] + '")']
    test('Compare' + i, octmat, args : octmat_args,
      suite : i, timeout: 60, workdir : meson.source_root() / 'tests')
  endif

endforeach