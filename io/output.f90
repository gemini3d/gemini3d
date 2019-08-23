submodule (io) output
use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
implicit none
contains


module procedure create_outdir
! subroutine create_outdir(outdir,infile,indatsize,indatgrid,flagdneu,sourcedir,flagprecfile,precdir,flagE0file,E0dir)
!! CREATES OUTPUT DIRECTORY, MOVES CONFIG FILES THERE AND GENERATES A GRID OUTPUT FILE

integer(c_int) :: ierr


!> MAKE A COPY OF THE INPUT DATA IN THE OUTPUT DIRECTORY
if (mkdir(outdir//'/inputs') /= 0) error stop 'error creating output directory'

ierr = copyfile(infile, outdir//'/inputs/')
if (ierr /= 0) error stop 'error copying configuration .ini to output directory'
ierr = copyfile(indatsize, outdir//'/inputs/')
if (ierr /= 0) error stop 'error copying input data size file to output directory'
ierr = copyfile(indatgrid, outdir//'/inputs/')
if (ierr /= 0) error stop 'error copying input parameters to output directory'
ierr = copyfile(indatfile, outdir//'/inputs/')
if (ierr /= 0) error stop 'error copying input parameters to output directory'

!MAKE COPIES OF THE INPUT DATA, AS APPROPRIATE
if (.false.) then
  if (flagdneu/=0) then
    ierr = mkdir(outdir//'/inputs/neutral_inputs')
    ierr = copyfile(sourcedir//'/*', outdir//'/inputs/neutral_inputs/')
  end if
  if (ierr /= 0) error stop 'error copying neutral input parameters to output directory'

  if (flagprecfile/=0) then
    ierr = mkdir(outdir//'/inputs/prec_inputs')
    ierr = copyfile(precdir//'/*', outdir//'/inputs/prec_inputs/')
  end if
  if (ierr /= 0) error stop 'error copying input precipitation parameters to output directory'

  if (flagE0file/=0) then
    ierr = mkdir(outdir//'/inputs/Efield_inputs')
    ierr = copyfile(E0dir//'/*', outdir//'/inputs/Efield_inputs/')
  end if
  if (ierr /= 0) error stop 'error copying input energy parameters to output directory'
endif

!NOW STORE THE VERSIONS/COMMIT IDENTIFIER IN A FILE IN THE OUTPUT DIRECTORY
! this can break on POSIX due to copying files in endless loop, commended out - MH
! ierr = mkdir(outdir//'/inputs/source/')
!if (ierr /= 0) error stop 'error creating input source parameter output directory'
!call execute_command_line('cp -r ./* '//outdir//'/inputs/source/', exitstat=ierr)
!if (ierr /= 0) error stop 'error creating input source parameter output directory'

call gitlog(outdir//'/gitrev.log')

call compiler_log(outdir//'/compiler.log')

end procedure create_outdir


subroutine gitlog(logpath)
!! logs git branch, hash to file

character(*), intent(in) :: logpath
integer :: ierr

!> write branch
call execute_command_line('git rev-parse --abbrev-ref HEAD > '// logpath, cmdstat=ierr)
if(ierr /= 0) then
  write(stderr, *) 'failed to log Git branch'
  return
endif

!> write hash
call execute_command_line('git rev-parse --short HEAD >> '// logpath, cmdstat=ierr)
if(ierr /= 0) then
  write(stderr, *) 'failed to log Git hash'
  return
endif

!> write changed filenames
call execute_command_line('git status --porcelain >> '// logpath, cmdstat=ierr)
if(ierr /= 0) then
  write(stderr, *) 'failed to log Git filenames'
  return
endif

end subroutine gitlog


subroutine compiler_log(logpath)

character(*), intent(in) :: logpath
integer :: u, ierr

open(newunit=u, file=logpath, status='unknown', action='write', iostat=ierr)
if(ierr /= 0) return

write(u,'(A,/)') compiler_version()
write(u,'(A)') compiler_options()

close(u)

end subroutine compiler_log

end submodule output