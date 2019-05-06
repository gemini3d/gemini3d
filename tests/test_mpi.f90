use, intrinsic :: iso_fortran_env
use mpi
implicit none

integer :: mrank, msize, vlen, ierr
character(MPI_MAX_LIBRARY_VERSION_STRING) :: version  ! allocatable not ok

print *,compiler_version(), compiler_options()

call MPI_INIT(ierr)
if (ierr /= 0) error stop 'mpi init error'

call MPI_COMM_RANK(MPI_COMM_WORLD, mrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, msize, ierr)
call MPI_GET_LIBRARY_VERSION(version, vlen, ierr)

call MPI_FINALIZE(ierr)
if (ierr /= 0) error stop 'mpi finalize error'

print '(A,I3,A,I3,A)', 'Image ', mrank, ' / ', msize, ':',version

print '(A25,A7)','type','value'
print *,'mpi_real',mpi_real
print *,'mpi_real8',mpi_real8

end program


! intel:
! type  value
! mpi_real  1275069468
! mpi_real8  1275070505

