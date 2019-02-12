project('GEMINI', 'fortran',
  default_options : ['default_library=static'])

REALBITS = '-DREALBITS='+get_option('realbits')
GLOW = get_option('GLOW')
HDF = get_option('HDF')
mpi_root = get_option('MPI_ROOT')
lapack_root = get_option('LAPACK_ROOT')
scalapack_root = get_option('SCALAPACK_ROOT')
mumps_root = get_option('MUMPS_ROOT')

add_global_arguments(['-O3'], language : 'fortran')
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_global_arguments('-march=native', '-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function', language : 'fortran')
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

mpiret = fc.run('use mpi; end', dependencies : mpi).returncode()
if mpiret != 0
  error('MPI does not link properly with compiler')
endif

if lapack_root == ''
  lapack = dependency('lapack', required : true)
else
  message('searching for LAPACK in ' + lapack_root)
  lapack = fc.find_library('lapack', required : true, dirs : lapack_root)
endif

# meson >= 0.50
# mumps = dependency('MUMPS') module_path : 'cmake/Modules')
if mumps_root == ''
  mumpsinc = include_directories('/usr/include')
  mumpslib = [fc.find_library('dmumps', required : true)]
else
  message('searching for MUMPS in ' + mumps_root)
  mumpsinc = include_directories(mumps_root / 'include')
  mumpslib = [fc.find_library('dmumps', required : true, dirs : mumps_root / 'lib'),
              fc.find_library('mumps_common', required : true, dirs : mumps_root / 'lib'),
              fc.find_library('scalapack', required : true, dirs : scalapack_root / 'lib')]
endif

if HDF  # see HDF5 branch, not yet fully ready.
 hdf5 = dependency('hdf5')
endif

mpi_src = files('numerical/mpimod/mpimod.F90', 'numerical/mpimod/mpisend.f90', 'numerical/mpimod/mpirecv.f90', 'numerical/mpimod/mpihalo.f90')

calculus = library('calculus', mpi_src, 'numerical/constants/phys_consts.F90', 'numerical/calculus/calculus.f90', 'numerical/grid/grid.f90', 'numerical/interpolation/interpolation.f90', 'vendor/lapack95/gbsv.F90',
  dependencies : mpi,
  fortran_args: REALBITS)

io = library('io','io/io.F90', 'io/expanduser.f90',
  link_with: calculus,
  fortran_args: REALBITS)

msis = library('msis', 'vendor/msis00/msis00_gfortran.f',
  fortran_args : oldflags)

neutral = library('neutral', 'neutral/neutral.f90', 'temporal/timeutils.f90',
  dependencies : mpi,
  link_with : [calculus, io, msis])

# --- ionization
ionization_src = files('ionization/ionization.f90', 'ionization/boundary_conditions/precipBCs_mod.f90')
if GLOW
  ionization_src += 'ionization/glow_run.f90'
  glow_proj = subproject('cglow')
  cglow = glow_proj.get_variable('cglow')
  ionization_deps = cglow
else
  ionization_src += 'ionization/glow_dummy.f90'
  ionization_deps = []
endif

ionization = library('ionization', ionization_src,
  link_with: [ionization_deps, io, calculus, neutral],
  dependencies : mpi)

# -- gemini exe
gemini = executable('gemini.bin', 'gemini.f90', 'temporal/temporal.f90', 'numerical/potential/potential_comm_mumps.f90', 'numerical/potential/potential_mumps.F90', 'collisions/collisions.f90', 'numerical/potential/boundary_conditions/potentialBCs_mumps.f90', 'multifluid/multifluid.f90',
'sources/sources.f90', 'numerical/advection/advec_mpi.f90', 'numerical/diffusion/diffusion.F90',

  link_with : [io,calculus, ionization, neutral],
  dependencies : [mpi, mumpslib, lapack],
  include_directories : mumpsinc,
  fortran_args: REALBITS)

# -- magcalc exe

magcalc = executable('magcalc.bin', 'magcalc.f90',
  link_with : [calculus, neutral, io],
  dependencies : mpi)

# --- TESTS ---
# meson test --wrap="mpirun -np 4"

oct = find_program('octave')
test2d = meson.build_root() / 'test2d'
test3d = meson.build_root() / 'test3d'
zenodo2d = meson.source_root() / '../simulations/zenodo2d'
eval2d = '''--eval compare_all("''' + test2d + '","' + zenodo2d + '")'

test('Gemini2D', gemini,
  args : [meson.source_root() / 'initialize/test2d/config.ini', test2d],
  workdir : meson.source_root(),
  timeout : 90,
  is_parallel : false)

test('Gemini3D', gemini,
  args : [meson.source_root() / 'initialize/test3d/config.ini', test3d],
  workdir : meson.source_root(),
  timeout : 180,
  is_parallel : false)

if oct.found()
  test('Compare2D', oct, args : ['-q', eval2d])
endif
