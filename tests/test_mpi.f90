program mpiTypes
use, intrinsic :: iso_fortran_env
use mpi_f08

implicit none

integer :: mrank, msize, vlen
character(MPI_MAX_LIBRARY_VERSION_STRING) :: version  ! allocatable not ok

print *,compiler_version(), compiler_options()

call MPI_INIT()
call MPI_COMM_RANK(MPI_COMM_WORLD, mrank)
call MPI_COMM_SIZE(MPI_COMM_WORLD, msize)
call MPI_GET_LIBRARY_VERSION(version, vlen)

call MPI_FINALIZE()

print '(A,I3,A,I3,A)', 'Image ', mrank, ' / ', msize, ':',version

print '(A25,A7)','type','value'
print *,'mpi_real',mpi_real
print *,'mpi_real8',mpi_real8

end program


! intel:
! type  value
! mpi_real  1275069468
! mpi_real8  1275070505

