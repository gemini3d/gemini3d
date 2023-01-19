program test_mpi

use mpi_f08
use, intrinsic :: iso_fortran_env, only : compiler_version, stderr=>error_unit

implicit none

character(6) :: argv

integer :: mrank, msize, vlen, ierr, N
character(MPI_MAX_LIBRARY_VERSION_STRING) :: version
!! allocatable character for version does not work

call get_command_argument(1, argv, status=ierr)
if(ierr /= 0) error stop "please specify number of MPI images (for checking)"
read(argv,*) N

call MPI_INIT()
call MPI_COMM_RANK(MPI_COMM_WORLD, mrank)
call MPI_COMM_SIZE(MPI_COMM_WORLD, msize)
call MPI_GET_LIBRARY_VERSION(version, vlen)

call MPI_FINALIZE()

if (N /= msize) then
  write(stderr,*) "ERROR: MPI image count from mpiexec:", N, "doesn't match mpi_comm_size:",msize
  error stop
endif

print '(A,I3,A,I3)', 'Image ', mrank, ' / ', msize-1
print *, 'MPI library version: ', version(:vlen)

if(mrank == 0) then
  print '(/,A,/)',compiler_version()
  print '(A12,A15)','type','value'
  print '(A12,I15)','mpi_real', MPI_REAL
  print '(A12,I15)','mpi_real8',MPI_REAL8
endif

end program

