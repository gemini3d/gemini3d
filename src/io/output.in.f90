submodule (io) output

use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options, real32, real64

implicit none (type, external)

contains


module procedure create_outdir
!! CREATES OUTPUT DIRECTORY, MOVES CONFIG FILES THERE AND GENERATES A GRID OUTPUT FILE

integer :: u, realbits
character(:), allocatable :: input_nml, output_dir, &
  git_version, branch, rev, porcelain, compiler, compiler_flags, exe
character(256) :: buf
integer :: mcadence
integer :: lid, lid2, lid3

namelist /files/ input_nml, output_dir,realbits
namelist /git/ git_version, branch, rev, porcelain
namelist /system/ compiler, compiler_flags, exe
namelist /grid/ lid, lid2, lid3, lx1, lx2, lx3, lx2all, lx3all
namelist /milestone/ mcadence

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
git_version = "@git_version@"
branch = "@git_branch@"
rev = "@git_rev@"
porcelain = "@git_porcelain@"

!> system namelist
compiler = trim(compiler_version())
compiler_flags = trim(compiler_options())
call get_command_argument(0, buf)
exe = trim(buf)

!> milestone namelist
mcadence=cfg%mcadence

lid = mpi_cfg%lid
lid2 = mpi_cfg%lid2
lid3 = mpi_cfg%lid3

!> let this crash the program if it can't write as an early indicator of output directory problem.
open(newunit=u, file=cfg%outdir // '/output.nml', status='unknown', action='write')
  write(u, nml=files)
  write(u, nml=git)
  write(u, nml=system)
  write(u, nml=grid)
  write(u, nml=milestone)
close(u)

end procedure create_outdir


subroutine compiler_log(logpath)

character(*), intent(in) :: logpath
integer :: u, ierr

open(newunit=u, file=logpath, status='unknown', action='write', iostat=ierr)
if(ierr /= 0) return

write(u,'(A,/,A)') compiler_version(), compiler_options()

close(u)

end subroutine compiler_log


end submodule output
