submodule (io) output
use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
implicit none
contains


module procedure create_outdir
! subroutine create_outdir(outdir,infile,indatsize,indatgrid,flagdneu,sourcedir,flagprecfile,precdir,flagE0file,E0dir)
!! CREATES OUTPUT DIRECTORY, MOVES CONFIG FILES THERE AND GENERATES A GRID OUTPUT FILE

integer(c_int) :: ierr


!> MAKE A COPY OF THE INPUT DATA IN THE OUTPUT DIRECTORY
if ( mkdir(outdir//'/inputs') /= 0 ) error stop 'error creating output directory'

if ( copyfile(infile, outdir//'/inputs/')  /= 0) error stop 'error copying configuration file to output directory'
if ( copyfile(indatsize, outdir//'/inputs/') /= 0) error stop 'error copying input data size file to output directory'
if ( copyfile(indatgrid, outdir//'/inputs/') /= 0) error stop 'error copying input grid to output directory'
if ( copyfile(indatfile, outdir//'/inputs/') /= 0) error stop 'error copying input data to output directory'

!MAKE COPIES OF THE INPUT DATA, AS APPROPRIATE
if (.false.) then
  if (flagdneu/=0) then
    ierr = mkdir(outdir//'/inputs/neutral_inputs')
    if ( copyfile(sourcedir//'/*', outdir//'/inputs/neutral_inputs/')  /= 0) error stop 'copy: neutral input => output dir'
  end if

  if (flagprecfile/=0) then
    ierr = mkdir(outdir//'/inputs/prec_inputs')
    if ( copyfile(precdir//'/*', outdir//'/inputs/prec_inputs/') /= 0) error stop 'copy: input precipitation => output dir'
  end if

  if (flagE0file/=0) then
    ierr = mkdir(outdir//'/inputs/Efield_inputs')
    if ( copyfile(E0dir//'/*', outdir//'/inputs/Efield_inputs/') /= 0) error stop 'copy input energy => output dir'
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