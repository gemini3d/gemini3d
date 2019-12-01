use mpi, only: mpi_init, mpi_finalize, mpi_comm_world
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit, i64=>int64, compiler_version, compiler_options

implicit none

#if REALBITS==32
include 'smumps_struc.h'
type (SMUMPS_STRUC) :: mumps_par
#else
include 'dmumps_struc.h'
type (DMUMPS_STRUC) :: mumps_par
#endif

integer :: ierr
integer(i64) :: i8

print *,compiler_version()

call mpi_init(ierr)
if (ierr /= 0) error stop 'mpi init error'
! Define a communicator for the package.
mumps_par%COMM = MPI_COMM_WORLD
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call simple_test(mumps_par)

call mpi_finalize(ierr)
if (ierr /= 0) error stop 'mpi finalize error'

contains


subroutine simple_test(mumps_par)

#if REALBITS==32
type (SMUMPS_STRUC), intent(inout) :: mumps_par
#else
type (DMUMPS_STRUC), intent(inout) :: mumps_par
#endif

call mumps_run(mumps_par)

!>  Define problem on the host (processor 0)
call read_input(mumps_par)

!>  Call package for solution
mumps_par%JOB = 6
call mumps_run(mumps_par)

!>  Solution has been assembled on the host
IF ( mumps_par%MYID .eq. 0 ) THEN
  print *, ' Solution is '
  print '(5F7.3)', mumps_par%RHS
END IF

!>  Deallocate user data
IF ( mumps_par%MYID .eq. 0 )THEN
  DEALLOCATE( mumps_par%IRN )
  DEALLOCATE( mumps_par%JCN )
  DEALLOCATE( mumps_par%A   )
  DEALLOCATE( mumps_par%RHS )
END IF

!>  Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
call mumps_run(mumps_par)

end subroutine simple_test


subroutine read_input(mumps_par)

#if REALBITS==32
type (SMUMPS_STRUC), intent(inout) :: mumps_par
#else
type (DMUMPS_STRUC), intent(inout) :: mumps_par
#endif

integer :: i, u


IF ( mumps_par%MYID == 0 ) THEN

open(newunit=u, file='input_simpletest_real', form='formatted', status='old', action='read')
READ(u,*) mumps_par%N
READ(u,*) mumps_par%NZ  ! NNZ for MUMPS > 5.1.0
ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
ALLOCATE( mumps_par%A( mumps_par%NZ ) )
ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
DO I8 = 1, mumps_par%NZ
  READ(u,*) mumps_par%IRN(I8),mumps_par%JCN(I8), mumps_par%A(I8)
END DO
DO I = 1, mumps_par%N
  READ(u,*) mumps_par%RHS(I)
END DO
close(u)

END IF

end subroutine read_input


subroutine mumps_run(mumps_par)

#if REALBITS==32
type (SMUMPS_STRUC), intent(inout) :: mumps_par
CALL SMUMPS(mumps_par)
#else
type (DMUMPS_STRUC), intent(inout) :: mumps_par
CALL DMUMPS(mumps_par)
#endif



IF (mumps_par%INFOG(1) < 0) THEN
  WRITE(stderr,'(A,A,I6,A,I9)') " ERROR: ", &
  "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),  "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)

  error stop
END IF

end subroutine mumps_run

end program
