project('GEMINI', 'fortran',
  default_options : ['default_library=static'])
  
REALBITS = '-DREALBITS='+get_option('realbits')
GLOW = get_option('GLOW')
HDF = get_option('HDF')

fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_global_arguments('-O3','-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function', language : 'fortran')
endif

mpi = dependency('mpi', language : 'fortran', required : true)
#mumps = dependency('mumps')#s, module_path : 'cmake/Modules')
mumps = fc.find_library('dmumps', required : true)
lapack = dependency('lapack', required : true)

if HDF  # see HDF5 branch, not yet fully ready.
 hdf5 = dependency('hdf5')
endif


calculus = library('calculus', 'numerical/constants/phys_consts.F90', 'numerical/calculus/calculus.f90', 'numerical/grid/grid.f90', 'numerical/mpimod/mpimod.F90', 'numerical/interpolation/interpolation.f90', 'vendor/lapack95/gbsv.F90',
  dependencies : mpi,
  fortran_args: REALBITS)

io = library('io','io/io.F90', 'io/expanduser.f90',
  link_with: calculus,
  fortran_args: REALBITS)
  
msis = library('msis', 'vendor/msis00/msis00_gfortran.f',
  fortran_args : '-w')

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
  dependencies : [mpi, mumps, lapack],
  include_directories: '/usr/include',
  fortran_args: REALBITS)
  
# -- magcalc exe

magcalc = executable('magcalc.bin', 'magcalc.f90',
  link_with : [calculus, neutral, io],
  dependencies : mpi)
  
# --- TESTS ---
# meson test --wrap="mpirun -np 4"

test('gemini 2D', gemini,
  args : [meson.source_root() / 'initialize/test2d/config.ini', meson.build_root() / 'test2d'],
  workdir : meson.source_root(),
  timeout : 240,
  is_parallel : false)
  
test('gemini 3D', gemini,
  args : [meson.source_root() / 'initialize/test3d/config.ini', meson.build_root() / 'test2d'],
  workdir : meson.source_root(),
  timeout : 600,
  is_parallel : false)
