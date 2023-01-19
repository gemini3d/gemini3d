program test_mumps

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit, compiler_version, compiler_options
use mpi_f08, only : mpi_init, mpi_comm_world,mpi_finalize
use mumps_interface, only: mumps_struc, mumps_exec

implicit none (type, external)

type(mumps_struc) :: mumps_par

call mpi_init()
! Define a communicator for the package.
mumps_par%COMM = MPI_COMM_WORLD%mpi_val
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call simple_test(mumps_par)

call mpi_finalize()

contains


subroutine simple_test(mumps_par)

type (MUMPS_STRUC), intent(inout) :: mumps_par

call mumps_run(mumps_par)

IF (mumps_par%MYID == 0) print *,compiler_version()

!>  Define problem on the host (processor 0)
IF (mumps_par%MYID == 0) call read_input(mumps_par)

!>  Call package for solution
mumps_par%JOB = 6
call mumps_run(mumps_par)

!>  Solution has been assembled on the host
IF (mumps_par%MYID == 0) print '(A,5F7.3)', ' Solution: ', mumps_par%RHS

!>  Deallocate user data
IF (mumps_par%MYID == 0) DEALLOCATE( mumps_par%IRN, mumps_par%JCN, mumps_par%A, mumps_par%RHS)

!>  Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
call mumps_run(mumps_par)

end subroutine simple_test


subroutine read_input(mumps_par)

type (MUMPS_STRUC), intent(inout) :: mumps_par

integer :: i, u
character(2048) :: argv
character(:), allocatable :: filename

integer :: N, NNZ
integer, allocatable :: IRN(:), JCN(:)
real, allocatable :: A(:), RHS(:)

namelist /shape/ N, NNZ
namelist /data/ IRN,JCN,A,RHS


call get_command_argument(1, argv, status=i)
if (i/=0) then
  filename = 'input_simpletest_real.nml'
else
  filename = trim(argv)
endif


open(newunit=u, file=filename, form='formatted', status='old', action='read')

read(u, nml=shape)
mumps_par%N = N
mumps_par%NNZ = NNZ
ALLOCATE(mumps_par%IRN(NNZ), IRN(NNZ), mumps_par%JCN(NNZ), JCN(NNZ), mumps_par%A(NNZ), A(NNZ), mumps_par%RHS(N), RHS(N))

read(u, nml=data)
mumps_par%IRN = IRN
mumps_par%JCN = JCN
mumps_par%A = A
mumps_par%RHS = RHS
close(u)

end subroutine read_input


subroutine mumps_run(mumps_par)

type (MUMPS_STRUC), intent(inout) :: mumps_par

call mumps_exec(mumps_par)

IF (mumps_par%INFOG(1) < 0) THEN
  WRITE(stderr,'(A,A,I6,A,I9)') " ERROR: ", &
  "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),  "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)

  error stop
END IF

end subroutine mumps_run

end program
