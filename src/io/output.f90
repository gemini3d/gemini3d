submodule (io) output

use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options, real32, real64

implicit none

contains


module procedure create_outdir
!! CREATES OUTPUT DIRECTORY, MOVES CONFIG FILES THERE AND GENERATES A GRID OUTPUT FILE

integer :: ierr, u, realbits
logical :: porcelain
character(:), allocatable :: branch, rev
character(256) :: buf

namelist /files/ infile,outdir,realbits
namelist /git/ branch, rev, porcelain

!> MAKE A COPY OF THE INPUT DATA IN THE OUTPUT DIRECTORY
ierr = mkdir(outdir//'/inputs')

ierr = copyfile(infile, outdir//'/inputs/')
ierr = copyfile(cfg%indatsize, outdir//'/inputs/')
ierr = copyfile(cfg%indatgrid, outdir//'/inputs/')
ierr = copyfile(cfg%indatfile, outdir//'/inputs/')

!MAKE COPIES OF THE INPUT DATA, AS APPROPRIATE
if (.false.) then
  if (cfg%flagdneu/=0) then
    ierr = mkdir(outdir//'/inputs/neutral_inputs')
    ierr = copyfile(cfg%sourcedir//'/*', outdir//'/inputs/neutral_inputs/')
  end if

  if (cfg%flagprecfile/=0) then
    ierr = mkdir(outdir//'/inputs/prec_inputs')
    ierr = copyfile(cfg%precdir//'/*', outdir//'/inputs/prec_inputs/')
  end if

  if (cfg%flagE0file/=0) then
    ierr = mkdir(outdir//'/inputs/Efield_inputs')
    ierr = copyfile(cfg%E0dir//'/*', outdir//'/inputs/Efield_inputs/')
  end if
endif

call gitlog(outdir // '/gitrev.log')

call compiler_log(outdir // '/compiler.log')

!! Log to output.nml
select case (wp)
case (real64)
  realbits = 64
case (real32)
  realbits = 32
case default
  error stop 'unknown real precision'
end select

branch = ''
rev = ''
porcelain = .false.
open(newunit=u, file=outdir// '/gitrev.log', status='old', action='read', iostat=ierr)
if(ierr==0) then
  read(u, '(A256)', iostat=ierr) buf
  if(ierr==0) branch = trim(buf)
  read(u, '(A256)', iostat=ierr) buf
  if(ierr==0) rev = trim(buf)
  read(u, '(A256)', iostat=ierr) buf
  if (len_trim(buf)==0 .or. is_iostat_end(ierr)) porcelain=.true.
  close(u)
endif

!> let this crash the program if it can't write as an early indicator of output directory problem.
open(newunit=u, file=outdir // '/output.nml', status='unknown', action='write')
  write(u, nml=files)
  write(u, nml=git)
close(u)

end procedure create_outdir


subroutine gitlog(logpath)
!! logs git branch, hash to file

character(*), intent(in) :: logpath
integer :: ierr

!> git branch --show-current requires Git >= 2.22, June 2019
call execute_command_line('git rev-parse --abbrev-ref HEAD > '// logpath, cmdstat=ierr)
if(ierr /= 0) then
  write(stderr, *) 'ERROR: failed to log Git branch'
  return
endif

!> write hash
call execute_command_line('git rev-parse --short HEAD >> '// logpath, cmdstat=ierr)
if(ierr /= 0) then
  write(stderr, *) 'ERROR: failed to log Git hash'
  return
endif

!> write changed filenames
call execute_command_line('git status --porcelain >> '// logpath, cmdstat=ierr)
if(ierr /= 0) then
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