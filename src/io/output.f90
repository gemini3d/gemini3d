submodule (io) output

use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options, stderr=>error_unit, real32, real64

implicit none

contains


module procedure create_outdir
! subroutine create_outdir(outdir,infile,indatsize,indatgrid,indatfile,flagdneu,sourcedir,flagprecfile,precdir,flagE0file,E0dir)
!! CREATES OUTPUT DIRECTORY, MOVES CONFIG FILES THERE AND GENERATES A GRID OUTPUT FILE

integer :: ierr


!> MAKE A COPY OF THE INPUT DATA IN THE OUTPUT DIRECTORY
ierr = mkdir(outdir//'/inputs')

ierr = copyfile(infile, outdir//'/inputs/')
ierr = copyfile(indatsize, outdir//'/inputs/')
ierr = copyfile(indatgrid, outdir//'/inputs/')
ierr = copyfile(indatfile, outdir//'/inputs/')

!MAKE COPIES OF THE INPUT DATA, AS APPROPRIATE
if (.false.) then
  if (flagdneu/=0) then
    ierr = mkdir(outdir//'/inputs/neutral_inputs')
    ierr = copyfile(sourcedir//'/*', outdir//'/inputs/neutral_inputs/')
  end if

  if (flagprecfile/=0) then
    ierr = mkdir(outdir//'/inputs/prec_inputs')
    ierr = copyfile(precdir//'/*', outdir//'/inputs/prec_inputs/')
  end if

  if (flagE0file/=0) then
    ierr = mkdir(outdir//'/inputs/Efield_inputs')
    ierr = copyfile(E0dir//'/*', outdir//'/inputs/Efield_inputs/')
  end if
endif

call gitlog(outdir // '/gitrev.log')

call compiler_log(outdir // '/compiler.log')

call realbits_log(outdir // '/realbits.log')

end procedure create_outdir


subroutine gitlog(logpath)
!! logs git branch, hash to file

character(*), intent(in) :: logpath
integer :: ierr

!> write branch
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


subroutine realbits_log(logpath)

character(*), intent(in) :: logpath
integer :: u, ierr

open(newunit=u, file=logpath, status='unknown', action='write', iostat=ierr)
if(ierr /= 0) return

select case (wp)
  case (real64)
    write(u,'(A)') '64'
  case (real32)
    write(u,'(A)') '32'
  case default
    write(u,'(A)') 'unknown'
end select

close(u)

end subroutine realbits_log

end submodule output