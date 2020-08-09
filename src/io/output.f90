submodule (io) output

use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options, real32, real64

implicit none (type, external)

contains


module procedure create_outdir
!! CREATES OUTPUT DIRECTORY, MOVES CONFIG FILES THERE AND GENERATES A GRID OUTPUT FILE

integer :: ierr, u, realbits
logical :: porcelain, exists
character(:), allocatable :: input_nml, output_dir, &
  git_version, branch, rev, compiler, compiler_flags, exe
character(256) :: buf
integer :: mcadence

namelist /files/ input_nml, output_dir,realbits
namelist /git/ git_version, branch, rev, porcelain
namelist /system/ compiler, compiler_flags, exe
namelist /grid/ lid, lid2, lid3, lx1, lx2, lx3, lx2all, lx3all
namelist /milestone/ mcadence

call gitlog(cfg%outdir // '/gitrev.log')

!> Log to output.nml

!> files namelist
input_nml = cfg%infile
output_dir = cfg%outdir
!! character namelist variables can't be assumed length, but can be allocatable.

select case (wp)
case (real64)
  realbits = 64
case (real32)
  realbits = 32
case default
  error stop 'unknown real precision'
end select

!> git namelist
git_version = ''
branch = ''
rev = ''
porcelain = .false.
open(newunit=u, file=cfg%outdir// '/gitrev.log', status='old', action='read', iostat=ierr)
if(ierr==0) then
  read(u, '(A256)', iostat=ierr) buf
  if(ierr==0) git_version = trim(buf)
  read(u, '(A256)', iostat=ierr) buf
  if(ierr==0) branch = trim(buf)
  read(u, '(A256)', iostat=ierr) buf
  if(ierr==0) rev = trim(buf)
  read(u, '(A256)', iostat=ierr) buf
  if (len_trim(buf)==0 .or. is_iostat_end(ierr)) porcelain=.true.
  close(u)
endif

!> system namelist
compiler = trim(compiler_version())
compiler_flags = trim(compiler_options())
call get_command_argument(0, buf)
exe = trim(buf)

!> milestone namelist
mcadence=cfg%mcadence

!> let this crash the program if it can't write as an early indicator of output directory problem.
open(newunit=u, file=cfg%outdir // '/output.nml', status='unknown', action='write')
  write(u, nml=files)
  write(u, nml=git)
  write(u, nml=system)
  write(u, nml=grid)
  write(u, nml=milestone)
close(u)

end procedure create_outdir


subroutine gitlog(logpath)
!! logs git branch, hash to file
!! this works on Windows and Unix

character(*), intent(in) :: logpath
integer :: i,j

call execute_command_line('git --version > ' // logpath, cmdstat=i, exitstat=j)
if(i /= 0 .or. j /= 0) then
  write(stderr, *) 'ERROR: Git not available'
  return
endif

!> git branch --show-current requires Git >= 2.22, June 2019
call execute_command_line('git rev-parse --abbrev-ref HEAD >> '// logpath, cmdstat=i, exitstat=j)
if(i /= 0 .or. j /= 0) then
  write(stderr, *) 'ERROR: failed to log Git branch'
  return
endif

!> write hash
call execute_command_line('git rev-parse --short HEAD >> '// logpath, cmdstat=i, exitstat=j)
if(i /= 0 .or. j /= 0) then
  write(stderr, *) 'ERROR: failed to log Git hash'
  return
endif

!> write changed filenames
call execute_command_line('git status --porcelain >> '// logpath, cmdstat=i, exitstat=j)
if(i /= 0 .or. j /= 0) then
  write(stderr, *) 'ERROR: failed to log Git filenames'
  return
endif

end subroutine gitlog


subroutine compiler_log(logpath)

character(*), intent(in) :: logpath
integer :: u, ierr

open(newunit=u, file=logpath, status='unknown', action='write', iostat=ierr)
if(ierr /= 0) return

write(u,'(A,/,A)') compiler_version(), compiler_options()

close(u)

end subroutine compiler_log


end submodule output
