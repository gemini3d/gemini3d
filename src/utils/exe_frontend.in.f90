module exe_frontend

use, intrinsic :: iso_c_binding, only : c_int
use, intrinsic :: iso_fortran_env, only : compiler_version, stderr=>error_unit, compiler_options
use phys_consts, only : wp
use gemini3d_config, only : gemini_cfg, read_configfile
use gemini3d_sysinfo, only : get_compiler_vendor
use filesystem, only : parent, assert_is_dir, expanduser, remove
use timeutils, only : date_filename,dateinc

implicit none (type, external)
private
public :: clean_output, cli_parser, get_Ncpu, help_gemini_bin, help_gemini_run, help_magcalc_bin, help_magcalc_run

interface !< cpu_count.cpp
integer(c_int) function cpu_count_c() bind(c, name="cpu_count")
import c_int
end function
end interface

contains


integer function cpu_count()
cpu_count = int(cpu_count_c())
end function


subroutine cli_parser(plan, Ncpu, path, exe, mpiexec, extra)

logical, intent(out) :: plan
integer, intent(out) :: Ncpu
character(:), allocatable, intent(out) :: path, exe, mpiexec, extra

character(1000) :: buf
integer :: argc, i, j, ierr, L

argc = command_argument_count()

call get_command_argument(1, buf, status=i)
if (i/=0) call help_run()

if(buf(1:1) == '-') then
!! not running sim, checking parameters
  select case(buf)
  case ('-h', '-help')
    call help_run()
  case ('-compiler')
    print '(A)', get_compiler_vendor()
  case ('-compiler_version')
    print '(A)', compiler_version()
  case ('-compiler_options')
    print '(A)', compiler_options()
  case ('-git')
    print '(A)', "@git_rev@"
  case ('-features')
    print '(A)', "@gemini_features@"
  case default
    write(stderr,*) "Gemini3D: unknown option: ", trim(buf)
    call help_run()
  end select

  stop
endif

!> simulation data directory
path = trim(buf)
call assert_is_dir(path)

plan = .false.
extra = ""
Ncpu = 0

do i = 2, argc
  call get_command_argument(i, buf)
  if(buf(1:1) /= "-") cycle  !< assume previous option was a -flag

  select case (buf)

  case ('-n')
    call get_command_argument(i+1, buf, length=L, status=ierr)
    if(ierr /= 0 .or. L==0 .or. buf(1:1) == "-") error stop trim(buf) // " -n missing parameter"
    read(buf, '(I6)') Ncpu
  case ('-exe')
    call get_command_argument(i+1, buf, length=L, status=ierr)
    if(ierr /= 0 .or. L==0 .or. buf(1:1) == "-") error stop trim(buf) // " -exe missing parameter"
    exe = find_exe(trim(buf))
  case ('-mpiexec')
    call get_command_argument(i+1, buf, length=L, status=ierr)
    if(ierr /= 0 .or. L==0 .or. buf(1:1) == "-") error stop "-mpiexec was specified without an executable path"
    mpiexec = find_mpiexec(trim(buf))
  case ('-plan')
    plan = .true.
  !> options passed to child executable
  case ('-dryrun', '-debug', '-nooutput')
    !! flags with no parameters
    extra = extra // ' ' // trim(buf)
  case ('-manual_grid')
    !! flags with two parameters
    extra = extra // ' ' // trim(buf)
    do j = 1,2
      call get_command_argument(i+j, buf, length=L, status=ierr)
      if(ierr /= 0 .or. L==0 .or. buf(1:1) == "-") error stop trim(buf) // " -manual_grid expected two parameters"
      extra = extra // ' ' // trim(buf)
    enddo
  case ('-start_time', '-end_time')
    !! flags with four parameters
    extra = extra // ' ' // trim(buf)
    do j = 1,4
      call get_command_argument(i+j, buf, length=L, status=ierr)
      if(ierr /= 0 .or. L==0 .or. buf(1:1) == "-") error stop trim(buf) // " -start_time expected four parameters"
      extra = extra // ' ' // trim(buf)
    enddo
  case default
    error stop "Gemini3D: unknown option: " // trim(buf)
  end select
end do

if(.not.allocated(exe)) exe = find_exe("")
if(.not.allocated(mpiexec)) mpiexec = find_mpiexec("")

end subroutine cli_parser


integer function get_Ncpu(N) result(Ncpu)

integer, intent(in) :: N
integer :: i
character(1000) :: buf

if (N > 0) then
  Ncpu = N
  print '(A,I0)', "run: CLI specified CPU count: ", Ncpu
  return
endif

call get_environment_variable("GEMINI_CPU", buf, status=i)
if (i==0) then
  read(buf,'(I6)', iostat=i) Ncpu
  if (i==0) print '(A,I0)', "run: GEMINI_CPU CPU job count: ", Ncpu
endif
if(i/=0) then
  call get_environment_variable("NSLOTS", buf, status=i)
  if (i==0) then
    read(buf,'(I6)', iostat=i) Ncpu
    if (i==0) print '(A,I0)', "run: SGE CPU job count: ", Ncpu
  endif
endif
if(i/=0) then
  call get_environment_variable("SLURM_NTASKS", buf, status=i)
  if(i==0) then
    read(buf, '(I6)', iostat=i) Ncpu
    if(i==0) print '(A,I0)', "run: SLURM CPU job count: ", Ncpu
  endif
endif
if (i/=0) then
  Ncpu = cpu_count()

  print '(A,I0)', 'gemini3d.run: detected CPU count: ', Ncpu
endif

end function get_Ncpu


function find_exe(name) result(exe)

character(*), intent(in) :: name

character(:), allocatable :: exe
character(8000) :: work, buf !< avoid quirks with reallocating, arbitrary lenght
logical :: exists

if(len_trim(name) > 0) then
  inquire(file=name, exist=exists)
  if(exists) then
    exe = name
    return
  endif
  inquire(file=name // '.exe', exist=exists)
  if(exists) then
    exe = name // '.exe'
    return
  endif
endif

call get_command_argument(0, buf)

if(len_trim(name) == 0) then
  if (index(buf, "gemini3d.run.debug") > 0) then
    work = "gemini.bin.debug"
  elseif (index(buf, "gemini3d.run") > 0) then
    work = "gemini.bin"
  elseif (index(buf, "magcalc.run.debug") > 0) then
    work = "magcalc.bin.debug"
  elseif (index(buf, "magcalc.run") > 0) then
    work = "magcalc.bin"
  else
    error stop "frontend: find_exe: please specify the program name you seek"
  endif

  work = parent(buf) // '/' // work

  inquire(file=work, exist=exists)
  if (exists) then
    exe = trim(work)
    return
  else
    inquire(file=trim(work) // '.exe', exist=exists)
    if(exists) then
      exe = trim(work) // '.exe'
      return
    endif
  endif
endif

error stop "gemini3d.run: did not find " // exe // " from " // name // " using " // trim(work) // &
  " : please specify path to MPI runnable executable with option 'gemini3d.run -exe path/to/my.bin'"

end function find_exe


function find_mpiexec(exe) result(mpiexec)

character(*), intent(in) :: exe
character(:), allocatable :: mpiexec

character(1000) :: buf
integer :: i, L

if(len_trim(exe) > 0) then
  if(check_mpiexec(expanduser(exe))) mpiexec = expanduser(exe)
else
  call get_environment_variable("MPI_ROOT", buf, length=L, status=i)
  if (i==0 .and. L>0) then
    if(check_mpiexec(expanduser(buf) // "/bin/mpiexec")) mpiexec = expanduser(buf) // "/bin/mpiexec"
  endif
endif

if(.not.allocated(mpiexec)) mpiexec = "mpiexec"

end function find_mpiexec


logical function check_mpiexec(exe) result(ok)

character(*), intent(in) :: exe

if(exe(1:1) == "-") error stop "gemini3d.run: -mpiexec was missing an executable, got: " // exe

inquire(file=exe, exist=ok)

if(ok) return

write(stderr,"(A,/,A)") "gemini3d.run: MPIexec file not found " // exe, &
"If simulation hangs or operates incorrectly, specify -mpiexec option or set MPI_ROOT environment variable."

end function check_mpiexec


subroutine help_gemini_bin() bind(C)

print '(/,A,/)', 'GEMINI-3D: gemini.bin ' // "@git_rev@"
print '(a,/,a)', 'by Matthew Zettergren', 'GLOW and auroral interfaces by Guy Grubbs'
print '(A)', 'build system and software engineering by Michael Hirsch'
print '(a,/,a)', 'Compiler vendor: '// get_compiler_vendor(), 'Compiler version: ' // compiler_version()
print '(/,A,/)', 'the first and only positional argument is simulation output directory.'
print '(A)', 'Optional arguments:'
print '(a,t25,a)', '-dryrun', 'allows quick check of first time step'
print '(a,t25,a)', '-manual_grid lx2 lx3', 'defines the number of MPI processes along x2 and x3.'
print '(t25,a)', '  If -manual_grid is not specified, the MPI processes are auto-assigned along x2 and x3.'
print '(a)', 'EOF: gemini.bin'

end subroutine help_gemini_bin


subroutine help_run()

character(1000) :: buf

call get_command_argument(0, buf)
if(index(buf, "gemini3d.run") > 0) call help_gemini_run()
if(index(buf, "magcalc.run") > 0) call help_magcalc_run()

error stop "help_run: unknown runner"

end subroutine help_run

subroutine help_gemini_run()

print '(/,a,/)', 'GEMINI-3D: gemini3d.run ' // "@git_rev@"
print '(a)', 'Compiler vendor: '// get_compiler_vendor()
print '(a)', 'Compiler version: ' // compiler_version()
print '(/,a,/)', 'the first and only positional argument is simulation output directory.'
print '(a)', 'Optional arguments:'
print '(a,t20,a)', '-plan', 'print MPI partition x2,x3 for given CPU count'
print '(a,t20,a)', '-dryrun', 'allows quick check of first time step'
print '(a,t20,a)', '-n', 'manually specify number of MPI images (default: auto-calculate MPI grid partitioning)'
print '(a,t20,a)', '-exe', 'specify path to gemini.bin (default: same directory as gemini3d.run)'
print '(a,t20,a)', '-mpiexec', 'specify path to mpiexec'
print '(a,t20,a)', '-compiler', 'tell the compiler family e.g. GNU, Intel, IntelLLVM.  This allows knowing what compiler was used.'
print '(a,t20,a)', '-compiler_version', 'like -compiler, and also tell the Fortran compiler version.'
print '(a,t20,a)', '-git', 'print git revision it was built from. This is not perfect, to be sure use fresh build directory.'
print '(a,t20,a)', '-compiler_options', 'print compiler flags used to build the executable.'
print '(a,t20,a)', '-features', 'print Gemini3D external features enabled as a string, machine/human-readable.'
stop 'EOF: gemini3d.run'

end subroutine help_gemini_run


subroutine help_magcalc_bin()

print '(/,A,/)', 'GEMINI-3D: magcalc.bin ' // "@git_rev@"
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'First two positional arguments must be input directory and fieldpoint file. Example:'
print '(/,A,/)', 'mpiexec -n 4 build/magcalc.bin test2d_fang test2d_fang/fieldpoint'
print '(A)', 'Optional arguments:'
print '(A)', '-dryrun option allows quick check of first time step'
stop 'EOF: magcalc.bin'

end subroutine help_magcalc_bin


subroutine help_magcalc_run()

print '(/,A,/)', 'GEMINI-3D: magcalc.run ' // "@git_rev@"
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A,/)', 'the first and only positional argument is simulation output directory.'
print '(A)', 'Optional arguments:'
print '(A)', '-plan  print MPI partition x2,x3 for given CPU count'
print '(A)', '-dryrun    allows quick check of first time step'
print '(A)', '-n   manually specify number of MPI images (default auto-calculate)'
print '(A)', '-exe   path to magcalc.bin'
print '(A)', '-mpiexec   path to mpiexec'
stop 'EOF: magcalc.run'

end subroutine help_magcalc_run


subroutine clean_output(path)

character(*), intent(in) :: path

type(gemini_cfg) :: cfg
integer, dimension(3) :: ymd
real(wp) :: UTsec
character(:), allocatable :: fn
logical :: exists

cfg%outdir = path
cfg%infile = path // '/inputs/config.nml'
inquire(file=cfg%infile, exist=exists)
if(.not.exists) error stop 'gemini3d.run: not a file: ' // cfg%infile

call read_configfile(cfg)

ymd = cfg%ymd0
UTsec = cfg%UTsec0

fn = date_filename(cfg%outdir, ymd, UTsec) // ".h5"

do
  !! new filename, add the 1 if it is the first
  fn = date_filename(cfg%outdir, ymd, UTsec) // ".h5"

  inquire(file=fn, exist=exists)
  if ( .not. exists ) exit
  !! last output file
  print *, 'delete: ', fn
  call remove(fn)

  !! next time
  call dateinc(cfg%dtout, ymd,UTsec)
end do

end subroutine clean_output


end module exe_frontend
