bits = get_option('realbits')
real_bits = '-DREALBITS=' + bits
arith = get_option('arith')
message('arith: ' + arith)

# === get compiler and OS

fc = meson.get_compiler('fortran')

os = host_machine.system()
if os == 'windows' and fc.get_id() != 'intel-cl'
  warning('On native Windows, only the Intel compilers (ifort) works at this time, because OpenMPI is not available for Windows. Consider Windows Subsystem for Linux.')
endif

# === compiler options
quiet = fc.first_supported_argument(['-w', '-nowarn', '/nowarn'])
impnone = fc.first_supported_argument(['-fimplicit-none', '-warn declarations', '/warn:declarations', '-Mdclchk'])
std18 = fc.first_supported_argument(['-std=f2018', '-stand=f18', '/stand:f18'])
bounds = fc.first_supported_argument(['-Werror=array-bounds', '-check bounds', '/check:bounds', '-Mbounds'])
if fc.get_id() == 'gcc'
  add_project_arguments('-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function',
    impnone, std18, language : 'fortran')
  if get_option('debug')
    add_project_arguments(bounds, '-fcheck=all', language: 'fortran')
  endif
elif fc.get_id() == 'intel'
  add_project_arguments(std18, impnone, language: 'fortran')
  if get_option('buildtype').contains('debug')
    # -heap-arrays didn't appear to be needed on Linux?
    add_project_arguments(bounds, language: 'fortran')
  endif
elif fc.get_id() == 'intel-cl'
  # /fpp is necessary for any Fortran preprocessing statement with Intel ifort on Windows
  # https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-fpp
   add_project_arguments('/fpp', std18, impnone, language: 'fortran')
   if get_option('buildtype').contains('debug')
    add_project_arguments(bounds, '/heap-arrays', language: 'fortran')
  endif
elif ['pgi', 'flang'].contains(fc.get_id())
  add_project_arguments(impnone, language: 'fortran')
  if get_option('buildtype').contains('debug')
    add_project_arguments(bounds, language: 'fortran')
  endif
endif

# -- MKL (in general, MKL can be used with intel, intel-cl, gcc or PGI)
mklroot = get_option('MKLROOT')

# == threads
threads = dependency('threads')

# == MPI
mpiexec = find_program('mpiexec')
mpi = dependency('mpi', language : 'fortran')
# this test needs to be .RUN() to catch some MPI configuration errors
if fc.run(files('../tests/test_mpi.f90'), dependencies : mpi, name: 'Fortran MPI').returncode() != 0
  error('mpi was not able to run a simple program')
endif

# --- Find Lapack

lapack_root = get_option('LAPACK_ROOT')
if lapack_root == ''
  blas = dependency('blas', required: false, disabler: true)
  if not blas.found()  # necessary for some systems incl. CentOS 7
    blas = fc.find_library('blas', required: false, disabler: true)
  endif
else
  blas = fc.find_library('blas', dirs: lapack_root, required: false, disabler: true)
endif
if blas.found() and not fc.links('x=dnrm2(1, [2.], 1); end', dependencies: blas, name: 'BLAS')
  blas = disabler()
endif

lapack = disabler()
if lapack_root == ''
  lapack = dependency('lapack', cmake_module_path: 'cmake/Modules', required: false, disabler: true)
  if not lapack.found()
    lapack = fc.find_library('lapack', required: false, disabler: true)
  endif
else
  lapack = fc.find_library('lapack', dirs: lapack_root, required: false, disabler: true)
endif
if lapack.found() and not fc.links('print *,disnan(0.); end', dependencies: lapack, name: 'LAPACK')
  lapack = disabler()
endif
if not lapack.found() or not blas.found()
  lapack_proj = subproject('lapack')
  lapack = lapack_proj.get_variable('lapack')
  blas = lapack_proj.get_variable('blas')
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
  if os == 'windows'
    blacsroot = mklroot/'lib/intel64_win'
  else
    blacsroot = mklroot/'lib/intel64'
  endif
  blacs += fc.find_library('mkl_blacs_intelmpi_lp64', dirs: blacsroot)
  if fc.get_id() == 'intel-cl'
    blacs += fc.find_library('impi', dirs: blacsroot)
    blacs += fc.find_library('iomp5md', dirs: blacsroot)
  else
    # blacs += fc.find_library('iomp5', dirs: blacsroot)  # link advisor cites, but not present on Linux...
  endif
endif


# --- Find Scalapack
# dependency('scalapack') not yet working in Meson--cmake->meson bug within Meson, yielding long meson-log.txt
# scalapack = dependency('scalapack', cmake_module_path : 'cmake/Modules')
scalapack = disabler()

if mklroot != ''
  scalapack_lib = fc.find_library('mkl_scalapack_lp64', dirs: blacsroot)
  scalapack = declare_dependency(dependencies: scalapack_lib, include_directories: mklroot/'include/intel64/lp64')
endif

if not scalapack.found()
  scalapack_root = get_option('SCALAPACK_ROOT')
  if scalapack_root != ''
    scalapack = fc.find_library('scalapack', dirs : scalapack_root / 'lib', required: false, disabler: true)
  else
    if os == 'linux'
      scalapack = fc.find_library('scalapack-openmpi', required: false, disabler: true)
      if not scalapack.found()
        scalapack = fc.find_library('scalapack', required: false, disabler: true)
      endif
    elif os == 'darwin'
      scalapack_root = '/usr/local/opt/scalapack'
      scalapack = fc.find_library('scalapack', dirs: scalapack_root / 'lib', required: false, disabler: true)
    endif
  endif
  if scalapack.found() and not fc.links('L = numroc(1,1,1,1,1); end', dependencies: scalapack, name: 'SCALAPACK')
    scalapack = disabler()
  endif
endif

if not scalapack.found()
  scalapack_proj = subproject('scalapack')
  scalapack = declare_dependency(link_with: scalapack_proj.get_variable('scalapack_lib'))
endif

# --- find M, 'bounds', umps
mumps = disabler()

if bits == '64'
  mumps_name = 'dmumps'
elif bits == '32'
  mumps_name = 'smumps'
else
  error('unknown precision ' + bits)
endif

mumps_root = get_option('MUMPS_ROOT')

if mumps_root == ''
  # mumps = dependency('MUMPS', cmake_module_path : 'cmake/Modules')
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

if mumps.found() and not fc.links(code, include_directories: include_directories(mumpsinc), dependencies: mumpslib, name: 'MUMPS')
  mumps = disabler()
endif
if not mumps.found() and scalapack.found()
  mumps_proj = subproject('mumps')
  mumps = mumps_proj.get_variable('mumps')
endif

# --- HDF5
# NOTE: Gemini HDF5 I/O is just a stub at the moment, undergoing development.
if get_option('hdf5').enabled()
  hdf5_proj = subproject('h5fortran')
  hdf5_interface = hdf5_proj.get_variable('hdf5_interface')
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

if run_command(python, '-c', 'import gemini').returncode() != 0
  warning('Please install PyGemini by "pip install -e ' + meson.source_root() + '"')
endif

#--- Octave

octave = find_program('octave', required: false, disabler: true)
# must do .found check since run_command(disabler()) has undefined returncode.
if octave.found() and run_command(octave, '--eval', 'assert(exist("validateattributes", "file")==2)').returncode() != 0
  message('Not using too-old Octave for comparison tests.')
  octave = disabler()
endif

#--- Matlab
# buggy for Matlab < R2019a, may hang "meson setup" silently for two minutes
matlab = find_program('matlab', required: get_option('matlab'), disabler: true)
