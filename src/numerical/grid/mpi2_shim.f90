module mpi2_shim
!! This module is a shim to workaround nuisance warnings in GCC >= 10 about interface mismatches

  use mpi, only : MPI_STATUS_SIZE, MPI_INTEGER

  implicit none (type, external)

  external :: mpi_send, mpi_recv

  private
  public :: mpi_send, mpi_send_int32_scalar, mpi_recv, mpi_recv_int32_scalar

contains

subroutine mpi_send_int32_scalar(BUF, DEST, TAG, COMM, IERROR)
  integer, intent(in) :: buf, dest, tag, comm, ierror
  call mpi_send(buf,1,MPI_INTEGER,dest,tag,COMM,ierror)
end subroutine mpi_send_int32_scalar

subroutine mpi_recv_int32_scalar(BUF, SOURCE, TAG, COMM, STATUS, IERROR)
  integer, intent(in) :: buf, source, tag, comm, STATUS(MPI_STATUS_SIZE), ierror
  call mpi_recv(buf,1,MPI_INTEGER,source,tag,COMM,ierror)
end subroutine mpi_recv_int32_scalar

end module mpi2_shim
