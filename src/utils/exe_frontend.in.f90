module exe_frontend

use, intrinsic :: iso_fortran_env, only : compiler_version, stderr=>error_unit
use config, only : get_compiler_vendor
use hwloc_ifc, only : get_cpu_count
use pathlib, only : parent, assert_directory_exists, expanduser

implicit none (type, external)

contains


subroutine cli_parser(plan, Ncpu, path, exe, mpiexec, extra)

logical, intent(out) :: plan
integer, intent(out) :: Ncpu
character(:), allocatable, intent(out) :: path, exe, mpiexec, extra

character(:), allocatable :: git_revision
character(1000) :: buf
integer :: argc, i

git_revision = "@git_rev@"

argc = command_argument_count()

call get_command_argument(1, buf, status=i)
if (i/=0) call help_run(git_revision)

if(buf(1:1) == '-') then
!! not running sim, checking parameters
  select case(buf)
  case ('-h', '-help')
    call help_run(git_revision)
  case ('-compiler')
    print '(A)', get_compiler_vendor()
  case ('-compiler_version')
    print '(A)', compiler_version()
  case ('-git')
    print '(A)', git_revision
  case default
    write(stderr,*) "Gemini3D: unknown option: ", trim(buf)
    call help_run(git_revision)
  end select

  stop
endif

!> simulation data directory
path = trim(buf)
call assert_directory_exists(path)

plan = .false.
exe = ""
mpiexec = ""
extra = ""
Ncpu = 0

do i = 2, argc
  call get_command_argument(i, buf)
  if(buf(1:1) /= "-") cycle  !< assume previous option was a -flag

  select case (buf)

  case ('-h', '-help')
    call help_run(git_revision)
  case ('-n')
    call get_command_argument(i+1, buf)
    read(buf, '(I6)') Ncpu
  case ('-exe', '-gemexe')
    !! FIXME: -gemexe is deprecated
    call get_command_argument(i+1, buf)
    exe = trim(buf)
  case ('-mpiexec')
    call get_command_argument(i+1, buf)
    mpiexec = trim(buf)
  case ('-dryrun')
    extra = '-dryrun'
  case ('-plan')
    plan = .true.
  case default
    write(stderr,*) "Gemini3D: unknown option: ", trim(buf)
  end select
end do

exe = find_exe(exe)
mpiexec = find_mpiexec(mpiexec)

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
  Ncpu = get_cpu_count()

  @apple_m1_workaround@

  print '(A,I0)', 'gemini3d.run: detected CPU count: ', Ncpu
endif

end function get_Ncpu


function find_exe(name, prog_name) result(exe)

character(*), intent(in) :: name
character(*), intent(in), optional :: prog_name

character(:), allocatable :: my, exe
logical :: exists
character(1000) :: buf

if(len_trim(name) > 0) then
  inquire(file=name, exist=exists)
  if(exists) then
    exe = name
    return
  endif
endif

call get_command_argument(0, buf)

if(present(prog_name)) then
  my = prog_name
else
  if(index(buf, "gemini3d.run") > 0) then
    my = "gemini.bin"
  elseif(index(buf, "magcalc.run") > 0) then
    my = "magcalc.bin"
  else
    error stop "please specify the program name you seek"
  endif
endif

if(len_trim(name) == 0) then
  if (len_trim(parent(buf)) > 0) then
    exe = trim(parent(buf)) // '/' // my
  else
    !! running from the same directory
    exe = my
  endif
endif
inquire(file=exe, exist=exists)
if(.not.exists) then
  inquire(file=exe // '.exe', exist=exists)
  if(exists) then
    exe = exe // '.exe'
  endif
endif

inquire(file=exe, exist=exists)
if(.not. exists) error stop "did not find " // exe // &
  " -- please specify path to MPI runnable executable with '-exe path/to/my.bin'"

end function find_exe


function find_mpiexec(exe) result(mpiexec)

character(*), intent(in) :: exe
character(:), allocatable :: mpiexec

character(1000) :: buf
integer :: i, L

mpiexec = ""

if(len_trim(exe) > 0) then
  if(check_mpiexec(expanduser(exe))) mpiexec = exe
else
  call get_environment_variable("MPI_ROOT", buf, length=L, status=i)
  if (i==0 .and. L>0) then
    if(check_mpiexec(expanduser(buf) // "/bin/mpiexec")) mpiexec = expanduser(buf) // "/bin/mpiexec"
  endif
endif

if(len_trim(mpiexec) == 0) mpiexec = "mpiexec"

end function find_mpiexec


logical function check_mpiexec(exe) result(ok)

character(*), intent(in) :: exe

inquire(file=exe, exist=ok)

if(ok) return

write(stderr,"(A)") "MPIexec file not found " // exe
write(stderr,"(A)") "If simulation hangs or operates incorrectly, specify -mpiexec option or set MPI_ROOT environment variable."

end function check_mpiexec


subroutine help_gemini_bin(git_revision)

character(*), intent(in) :: git_revision

print '(/,A,/)', 'GEMINI-3D: gemini.bin ' // git_revision
print '(A)', 'by Matthew Zettergren'
print '(A)', 'GLOW and auroral interfaces by Guy Grubbs'
print '(A)', 'build system and software engineering by Michael Hirsch'
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'must specify simulation output directory. Example:'
print '(/,A,/)', '  mpiexec -np 4 build/gemini.bin /path/to/simulation_outputs'
print '(A)', '-dryrun    allows quick check of first time step'
print '(A)', '-manual_grid lx2 lx3    defines the number of MPI processes along x2 and x3.'
print '(A)', '  If -manual_grid is not specified, the MPI processes are auto-assigned along x2 and x3.'
stop 'EOF: gemini.bin'

end subroutine help_gemini_bin


subroutine help_run(git_revision)

character(*), intent(in) :: git_revision
character(1000) :: buf

call get_command_argument(0, buf)
if(index(buf, "gemini3d.run") > 0) call help_gemini_run(git_revision)
if(index(buf, "magcalc.run") > 0) call help_magcalc_run(git_revision)

error stop "help_run: unknown runner"

end subroutine help_run

subroutine help_gemini_run(git_revision)

character(*), intent(in) :: git_revision

print '(/,A,/)', 'GEMINI-3D: gemini3d.run ' // git_revision
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'must specify simulation output directory. Example:'
print '(/,A,/)', '  build/gemini3d.run /path/to/simulation_outputs'
print '(A)', '-plan  print MPI partition x2,x3 for given CPU count'
print '(A)', '-dryrun    allows quick check of first time step'
print '(A)', '-n   manually specify number of MPI images (default auto-calculate)'
print '(A)', '-exe   path to gemini.bin'
print '(A)', '-mpiexec   path to mpiexec'
stop 'EOF: gemini3d.run'

end subroutine help_gemini_run


subroutine help_magcalc_bin(git_revision)
character(*), intent(in) :: git_revision

print '(/,A,/)', 'GEMINI-3D: magcalc.bin ' // git_revision
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'must specify input directory and fieldpoint file. Example:'
print '(/,A,/)', 'mpiexec -n 4 build/magcalc.bin test2d_fang test2d_fang/fieldpoint'
print '(A)', '-dryrun option allows quick check of first time step'
stop 'EOF: magcalc.bin'
end subroutine help_magcalc_bin

subroutine help_magcalc_run(git_revision)

character(*), intent(in) :: git_revision

print '(/,A,/)', 'GEMINI-3D: magcalc.run ' // git_revision
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'must specify simulation output directory. Example:'
print '(/,A,/)', '  build/magcalc.run /path/to/simulation_outputs'
print '(A)', '-plan  print MPI partition x2,x3 for given CPU count'
print '(A)', '-dryrun    allows quick check of first time step'
print '(A)', '-n   manually specify number of MPI images (default auto-calculate)'
print '(A)', '-exe   path to magcalc.bin'
print '(A)', '-mpiexec   path to mpiexec'
stop 'EOF: magcalc.run'

end subroutine help_magcalc_run

end module exe_frontend
