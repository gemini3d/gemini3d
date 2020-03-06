use, intrinsic :: iso_fortran_env
use mpi
implicit none

integer :: mrank, msize, vlen, ierr
character(MPI_MAX_LIBRARY_VERSION_STRING) :: version
!! allocatable character for version does not work

print '(/,A,/)',compiler_version()

call MPI_INIT(ierr)
if (ierr /= 0) error stop 'mpi init error'

call MPI_COMM_RANK(MPI_COMM_WORLD, mrank, ierr)
if (ierr /= 0) error stop 'mpi comm_rank error'
call MPI_COMM_SIZE(MPI_COMM_WORLD, msize, ierr)
if (ierr /= 0) error stop 'mpi comm_size error'
call MPI_GET_LIBRARY_VERSION(version, vlen, ierr)
if (ierr /= 0) error stop 'mpi get_library_version error'

call MPI_FINALIZE(ierr)
if (ierr /= 0) error stop 'mpi finalize error'

print '(A,I3,A,I3,A)', 'Image ', mrank, ' / ', msize, ': MPI library version:', trim(version)

print '(A12,A15)','type','value'
print '(A12,I15)','mpi_real',mpi_real
print '(A12,I15)','mpi_real8',mpi_real8

end program


! intel:
! type  value
! mpi_real  1275069468
! mpi_real8  1275070505
