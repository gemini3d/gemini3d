bits = get_option('realbits')
real_bits = '-DREALBITS=' + bits
arith = get_option('arith')

nomalloc_env = environment({'MALLOC_PERTURB_': '0'})

# === get compiler and OS

fc = meson.get_compiler('fortran')

os = host_machine.system()

# === compiler options
oldflags = fc.get_supported_arguments(['-std=legacy', '-w', '-fno-implicit-none', '-Wno-error=array-bounds', '-fcheck=no-all',
  '-nocheck', '/nocheck', '-nowarn', '/nowarn', '-nostand', '/nowarn'])
std18 = fc.first_supported_argument(['-std=f2018', '-stand=f18', '/stand:f18'])
if fc.get_id() == 'gcc'
  add_project_arguments('-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function',
    '-fimplicit-none', std18, language : 'fortran')
  if get_option('buildtype').contains('debug')
    add_project_arguments('-Werror=array-bounds', '-fcheck=all', language: 'fortran')
  endif
elif fc.get_id() == 'intel'
  add_project_arguments('-warn', '-diag-disable 5268', '-diag-disable 7712', std18, language: 'fortran')
  if get_option('buildtype').contains('debug')
    # -heap-arrays didn't appear to be needed on Linux?
    add_project_arguments('-check bounds', language: 'fortran')
  endif
elif fc.get_id() == 'intel-cl'
  # /fpp is necessary for any Fortran preprocessing statement with Intel ifort on Windows
  # https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-fpp
  add_project_arguments('/warn', '/Qdiag-disable:5268', '/Qdiag-disable:7712', '/fpp', std18, language: 'fortran')
  if get_option('buildtype').contains('debug')
    add_project_arguments('/check:bounds', '/heap-arrays', language: 'fortran')
  endif
elif ['pgi', 'flang'].contains(fc.get_id())
  if get_option('buildtype').contains('debug')
    add_project_arguments('-Mbounds', language: 'fortran')
  endif
endif

# -- MKL (in general, MKL can be used with intel, intel-cl, gcc or PGI)
mklroot = get_option('mklroot')
if ['intel', 'intel-cl'].contains(fc.get_id()) and mklroot == ''
  missing_opt = os == 'windows' ? '-Dmklroot=%MKLROOT%   or   -Dmklroot=$env:mklroot' : '-Dmklroot=$MKLROOT'
  error('Intel compiler needs MKL via MKLROOT environment variable. Please use option ' + missing_opt)
endif
mkllib_root = mklroot/'lib/intel64'
# os == 'windows' ? mklroot/'lib/intel64_win' : mklroot/'lib/intel64'

# == threads
threads = dependency('threads')

# == MPI
mpiexec = find_program('mpiexec')
mpi = dependency('mpi', language : 'fortran')
# this test needs to be .RUN() to catch some MPI configuration errors
if fc.run(files('../tests/test_mpi.f90'), dependencies : mpi, name: 'Fortran MPI').returncode() != 0
  error('mpi was not able to run a simple program')
endif

# == Intel OpenMP
intelMP_lib = []
if mklroot != ''
  n = os == 'windows' ? 'iomp5md' : 'iomp5'
  intelMP_lib = fc.find_library(n, dirs: mklroot / '../compiler/lib/intel64')
endif

# == Lapack
if mklroot == ''
  blas = dependency('blas', required: false, disabler: true)
  if not blas.found()  # necessary for some systems incl. CentOS 7
    blas = fc.find_library('blas', required: false, disabler: true)
  endif
  if not fc.links('x=dnrm2(1, [2.], 1); end', dependencies: blas, name: 'BLAS')
    blas = disabler()
  endif

  lapack = dependency('lapack', required: false, disabler: true)
  if not lapack.found()
    lapack = fc.find_library('lapack', required: false, disabler: true)
  endif
  if not fc.links('print *,disnan(0.); end', dependencies: lapack, name: 'LAPACK')
    lapack = disabler()
  endif
  if not lapack.found() or not blas.found()
    lapack_proj = subproject('lapack')
    lapack = declare_dependency(link_with: lapack_proj.get_variable('lapack'))
    blas = declare_dependency(link_with: lapack_proj.get_variable('blas'))
  endif
else
  # MKL
  lapack = []
  n1 = fc.get_id() == 'gcc' ? 'mkl_gf_lp64' : 'mkl_intel_lp64'
  names = [n1, 'mkl_intel_thread', 'mkl_core']
  foreach n : names
    lapack += fc.find_library(n, dirs: mkllib_root)
  endforeach
  lapack += intelMP_lib
  blas = lapack
endif

# --- intel MPI lib for Windows
intelmpi_lib = []
if mklroot != '' and os == 'windows'
  intelmpi_lib = fc.find_library('impi', dirs: mklroot / '../mpi/intel64/lib/release')
endif

# --- Find BLACS
# some systems need Blacs explicitly, but most don't, so don't pass disabler() if not found.
blacs = []
if mklroot == ''
  blacs = fc.find_library('blacs-openmpi', required: false)
  if not blacs.found()
    blacs = []
  endif
else
  blacs += fc.find_library('mkl_blacs_intelmpi_lp64', dirs: mkllib_root)
  blacs += intelmpi_lib
endif


# --- Find Scalapack
scalapack = dependency('scalapack', cmake_module_path: meson.source_root() /'cmake/Modules', required: false, disabler: true)

if not scalapack.found()
  foreach n : ['scalapack-openmpi', 'scalapack']
    scalapack = fc.find_library(n, required: false, disabler: true)
    if scalapack.found()
      break
    endif
  endforeach
endif

# FIXME: make real32 test
if fc.run(files('../tests/test_scalapack_d.f90'), dependencies: [scalapack, mpi], name: 'SCALAPACK').returncode() != 0
  # .links() is NOT an adequate test
  scalapack = disabler()
endif

if not scalapack.found()
  scalapack_proj = subproject('scalapack')
  scalapack = declare_dependency(link_with: scalapack_proj.get_variable('scalapack_lib'))
endif

# --- find Mumps
mumps = disabler()

if mklroot == ''
  if bits == '64'
    mumps_name = 'dmumps'
  elif bits == '32'
    mumps_name = 'smumps'
  else
    error('unknown precision ' + bits)
  endif

  mumps_root = get_option('MUMPS_ROOT')
  if mumps_root == ''
    if os == 'linux'
      mumps_root = '/usr'
      mumpslib = [fc.find_library(mumps_name, required : false, disabler: true),
                  fc.find_library('mumps_common', required: false, disabler: true)]
    elif os == 'darwin'
      # assuming homebrew
      mumps_root = '/usr/local/opt/mumps'
      mumpslib = [fc.find_library(mumps_name, dirs: mumps_root / 'lib', required : false, disabler: true),
                    fc.find_library('mumps_common', dirs: mumps_root / 'lib', required: false, disabler: true)]
    else
      mumpslib = [disabler(), disabler()]
    endif
  else
    mumpslib = [fc.find_library(mumps_name, dirs : mumps_root, required : false, disabler: true),
                fc.find_library('mumps_common', dirs : mumps_root, required : false, disabler: true)]
  endif

  mumpsinc = mumps_root / 'include'
  if mumpslib[0].found() and mumpslib[1].found()
    mumps = declare_dependency(include_directories: include_directories(mumpsinc), dependencies: mumpslib)
  endif

  code='''
  include 'dmumps_struc.h'
  type(DMUMPS_STRUC) :: mumps_par
  end'''
  # NOTE: leave this "if mumps.found()" so that mumpsinc is sure to be defined
  if mumps.found() and not fc.links(code, include_directories: include_directories(mumpsinc), dependencies: mumpslib, name: 'MUMPS')
    mumps = disabler()
  endif
endif

if not mumps.found()
  # MKL will always build MUMPS too
  mumps_proj = subproject('mumps')
  mumps = mumps_proj.get_variable('mumps')
endif

# --- NetCDF4
if get_option('netcdf').enabled()
  netcdf_proj = subproject('nc4fortran')
  netcdf_interface = netcdf_proj.get_variable('netcdf_interface')
  assert(netcdf_interface.found(), 'NetCDF I/O requested but NetCDF library not available')
else
  netcdf_interface = disabler()
endif

# --- HDF5
if not netcdf_interface.found() and get_option('hdf5').enabled()
  hdf5_proj = subproject('h5fortran')
  hdf5_interface = hdf5_proj.get_variable('hdf5_interface')
  assert(hdf5_interface.found(), 'HDF5 I/O requested but HDF5 library not available')
else
  hdf5_interface = disabler()
endif

# --- Python

pymod = import('python')
python = pymod.find_installation('python3')

message('Using Python ' + python.language_version())

if python.language_version().version_compare('< 3.6')
  error('Python ' + python.language_version() + ' is too old, you must have at least Python 3.6. Try "conda activate" or "module load python3" or similar.')
endif

# not needed at the moment
# is_ci = run_command(python, meson.source_root() / 'script_utils/meson_get_environment_variable.py', 'CI').stdout().to_lower() == 'true'

if run_command(python, '-c', 'import gemini3d').returncode() != 0
  message('Please install PyGemini by "pip install -e ' + meson.source_root() + '"')
endif

#--- Octave

octave = find_program('octave', required: false, disabler: true)
# must do .found check since run_command(disabler()) has undefined returncode.
if run_command(octave, '--eval', 'assert(exist("validateattributes", "file")==2)').returncode() != 0
  message('SKIP: Octave version is too old for comparison tests.')
  octave = disabler()
endif

#--- Matlab
# buggy for Matlab < R2019a, may hang "meson setup" silently for two minutes
matlab = find_program('matlab', required: get_option('matlab'), disabler: true)
