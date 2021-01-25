program get_mpi_count
!! for use from terminal/CMake, computes optimal MPI count
!! for a particular simulation

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use phys_consts, only : wp
use config, only : read_configfile, gemini_cfg
use autogrid, only : grid_auto, max_mpi
use reader, only: get_simsize3
use timeutils, only : date_filename,dateinc
use pathlib, only : get_suffix

implicit none (type, external)

integer :: i, lx1, lx2all, lx3all, lid, lid2, lid3, Ncpu, argc
character(1000) :: buf
character(:), allocatable :: path, gem_exe, cmd, mpiexec, extra
logical :: exists

argc = command_argument_count()

if(argc < 2) error stop 'specify top-level simulation path and cpu count'

call get_command_argument(1, buf, status=i)
path = trim(buf)

call get_command_argument(2, buf)
read(buf, '(I6)') Ncpu

call get_simsize3(path // '/inputs', lx1, lx2all, lx3all)

lid = max_mpi(lx2all, lx3all, Ncpu)

!> checks consistency
call grid_auto(lx2all, lx3all, lid, lid2, lid3)

if(argc < 3) then
  write(buf, '(I6,A)') lid, ' CPU cores would be used. Give path to gemini.bin to run simulation'
  stop trim(buf)
end if

call get_command_argument(3, buf)
gem_exe = trim(buf)
inquire(file=gem_exe, exist=exists)
if(.not.exists) error stop gem_exe // ' is not a file.'

if (argc >= 4) then
  call get_command_argument(4, buf)
  mpiexec = trim(buf)
  inquire(file=mpiexec, exist=exists)
  if(.not.exists) error stop mpiexec // ' is not a file.'
else
  mpiexec = 'mpiexec'
endif

if (argc >= 5) then
  call get_command_argument(5, buf)
  extra = trim(buf)
else
  extra = ''
endif

!> remove old output files
call clean_output(path)

!> run gemini.bin
write(buf, '(A1,A,A,I6,1X,A,1X,A,1X,A)') '"',mpiexec, '" -n ', lid, gem_exe, path, extra
cmd = trim(buf)
print *, cmd
call execute_command_line(cmd, exitstat=i)

if(i/=0) error stop 'gemini.bin run failure'


contains


subroutine clean_output(path)

character(*), intent(in) :: path

type(gemini_cfg) :: cfg
integer, dimension(3) :: ymd
real(wp) :: UTsec
character(:), allocatable :: fn, suffix
logical :: exists

cfg%outdir = path
cfg%infile = path // '/inputs/config.nml'
inquire(file=cfg%infile, exist=exists)
if(.not.exists) error stop 'mpi_cli: not a file: ' // cfg%infile

call read_configfile(cfg)

ymd = cfg%ymd0
UTsec = cfg%UTsec0

suffix = get_suffix(cfg%indatsize)
fn = date_filename(cfg%outdir, ymd, UTsec) // suffix

do
  !! new filename, add the 1 if it is the first
  fn = date_filename(cfg%outdir, ymd, UTsec) // suffix

  inquire(file=fn, exist=exists)
  if ( .not. exists ) exit
  !! last output file
  print *, 'delete: ', fn
  call unlink(fn)

  !! next time
  call dateinc(cfg%dtout, ymd,UTsec)
end do

end subroutine clean_output


subroutine unlink(path)
character(*), intent(in) :: path
integer :: i
logical :: e

inquire(file=path, exist=e)
if (.not.e) return

open(newunit=i, file=path, status='old')
close(i, status='delete')
end subroutine unlink


end program
