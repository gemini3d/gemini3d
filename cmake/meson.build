fc = meson.get_compiler('fortran')

os = host_machine.system()
if os == 'windows' and fc.get_id() != 'intel-cl'
  warning('On native Windows, only the Intel compilers (ifort) works at this time, because OpenMPI is not available for Windows. Consider Windows Subsystem for Linux.')
endif

# -- compiler options
quiet = fc.first_supported_argument(['-w', '-nowarn', '/nowarn'])
impnone = fc.first_supported_argument(['-fimplicit-none', '-warn declarations', '/warn:declarations', '-Mdclchk'])
std18 = fc.first_supported_argument(['-std=f2018', '-stand f18', '/stand:f18'])
bounds = fc.first_supported_argument(['-Werror=array-bounds', '-check bounds', '/check:bounds', '-Mbounds'])
if fc.get_id() == 'gcc'
  add_project_arguments('-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function',
    impnone, std18, language : 'fortran')
  if buildtype.contains('debug')
    add_project_arguments(bounds, '-fcheck=all', language: 'fortran')
  endif
elif fc.get_id() == 'intel'
  add_project_arguments(std18, impnone, language: 'fortran')
  if buildtype.contains('debug')
    # -heap-arrays didn't appear to be needed on Linux?
    add_project_arguments(bounds, language: 'fortran')
  endif
elif fc.get_id() == 'intel-cl'
  # /fpp is necessary for any Fortran preprocessing statement with Intel ifort on Windows
  # https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-fpp
   add_project_arguments('/fpp', std18, impnone, language: 'fortran')
   if buildtype.contains('debug')
    add_project_arguments(bounds, '/heap-arrays', language: 'fortran')
  endif
elif ['pgi', 'flang'].contains(fc.get_id())
  add_project_arguments(impnone, language: 'fortran')
  if buildtype.contains('debug')
    add_project_arguments(bounds, language: 'fortran')
  endif
endif

threads = dependency('threads')

# == MPI
mpiexec = find_program('mpiexec', required: false, disabler: true)
mpi = dependency('mpi', language : 'fortran')
if mpi.found() and not fc.links('use mpi; end', dependencies : mpi, name: 'Fortran MPI')
  mpi = disabler()
endif

# --- Find Lapack
# we don't use fallback because some systems have screwed-up libs factory installed
lapack = disabler()
blas = disabler()
system_blas = get_option('system_blas')
lapack_root = get_option('LAPACK_ROOT')
if system_blas
  blas = dependency('blas', required: false, disabler: true)
  if not blas.found()  # necessary for some systems incl. CentOS 7
    blas = fc.find_library('blas', dirs: lapack_root, required: false, disabler: true)
  endif
endif
if blas.found() and not fc.links('x=dnrm2(1, [2.], 1); end', dependencies: blas, name: 'BLAS')
  blas = disabler()
endif

lapack = disabler()
if system_blas and blas.found()
  lapack = dependency('lapack', cmake_module_path: 'cmake/Modules', required: false, disabler: true)
  if not lapack.found()
    lapack = fc.find_library('lapack', dirs: lapack_root, required: false, disabler: true)
  endif
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
blacs = fc.find_library('blacs-openmpi', required: false)
if not blacs.found()
  blacs = []
endif

# --- Find Scalapack
scalapack = disabler()

if system_blas
  # dependency('scalapack') not yet working in Meson--cmake->meson bug within Meson, yielding long meson-log.txt
  # scalapack = dependency('scalapack', cmake_module_path : 'cmake/Modules')
  scalapack_root = get_option('SCALAPACK_ROOT')
  if scalapack_root == ''  # system
    scalapack = fc.find_library('scalapack-openmpi', required: false, disabler: true)
    if not scalapack.found()
      scalapack = fc.find_library('scalapack', required: false, disabler: true)
    endif
  else
    message('searching for SCALAPACK in ' + scalapack_root)
    scalapack = fc.find_library('scalapack', dirs : scalapack_root / 'lib', required: false, disabler: true)
  endif
endif
if scalapack.found() and not fc.links('L = numroc(1,1,1,1,1); end', dependencies: scalapack, name: 'SCALAPACK')
  scalapack = disabler()
endif
if not scalapack.found() and mpiexec.found()
  scalapack_proj = subproject('scalapack')
  scalapack = declare_dependency(link_with: scalapack_proj.get_variable('scalapack_lib'))
endif

# --- find Mumps
mumps = disabler()

if bits == '64'
  mumps_name = 'dmumps'
elif bits == '32'
  mumps_name = 'smumps'
else
  error('unknown precision ' + bits)
endif
if system_blas
  # mumps = dependency('MUMPS', cmake_module_path : 'cmake/Modules')
  mumps_root = get_option('MUMPS_ROOT')
  if mumps_root == ''
    mumpsinc = []
    mumpslib = [fc.find_library(mumps_name, required : false, disabler: true),
                fc.find_library('mumps_common', required: false, disabler: true)]
    if os == 'linux'
      mumpsinc = '/usr/include'
    elif os == 'darwin'  # assuming homebrew
      mumpinc = '/usr/local/opt/mumps'
    endif
  else
    message('searching for MUMPS in ' + mumps_root)
    mumpsinc = mumps_root / 'include'
    mumpslib = [fc.find_library(mumps_name, dirs : mumps_root, required : false, disabler: true),
                fc.find_library('mumps_common', dirs : mumps_root, required : false, disabler: true)]
  endif
  if mumpslib[0].found() and mumpslib[1].found()  # has to be in if statement in case mumpsinc does not exist
    mumps = declare_dependency(include_directories: include_directories(mumpsinc), dependencies: mumpslib)
  endif
endif
code='''
include 'dmumps_struc.h'
type(DMUMPS_STRUC) :: mumps_par
end'''
if mumps.found() and not fc.links(code, include_directories: mumpsinc, dependencies: mumpslib, name: 'MUMPS')
  mumps = disabler()
endif
if not mumps.found() and scalapack.found()
  mumps_proj = subproject('mumps')
  mumps = mumps_proj.get_variable('mumps')
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

matlab = disabler()
if get_option('usematlab')  # buggy for Matlab < R2019a, may hang "meson setup" silently for two minutes
  matlab = find_program('matlab', required: false, disabler: true)
  # must do .found check since run_command(disabler()) has undefined returncode.
  if matlab.found() and run_command(python, meson.source_root() / 'script_utils/meson_matlab_has_batch.py').returncode() != 0
    matlab = disabler()
  endif
endif
