!------ ADAPTED FROM LAPACK_95 by Guy Grubbs --------
module lapack95
use, intrinsic :: iso_fortran_env, only: real32, real64
implicit none

#if REALBITS==32
integer, parameter :: wp=real32
#elif REALBITS==64
integer, parameter :: wp=real64
#else
error stop "realbits must be 32 or 64"
#endif

contains

subroutine gbsv(A,B,KL,IPIV,INFO)
  real(wp), dimension(:,:), intent(inout) :: A
  real(wp), dimension(:), intent(inout) :: B
  integer, intent(in), optional :: KL
  integer, dimension(:), intent(out), optional, target :: IPIV
  integer, intent(out), optional :: INFO
  integer :: LINFO, ISTAT, ISTAT1, SIPIV, LDA, N, NRHS, LKL, KU
  integer, dimension(:), pointer :: LPIV
  intrinsic size, present

  LINFO = 0; ISTAT = 0;
  LDA = size(A,1); N = size(A,2); NRHS = 1;

  if ( present(KL) ) then
    LKL = KL
  else
    LKL = (LDA-1)/3
  end if

  if ( present(IPIV) ) then
    SIPIV = size(IPIV)
  else
    SIPIV = N
  end if

  if ( present(IPIV) ) then
    LPIV => IPIV
  else
    allocate(LPIV(N))
  end if

  if ( ISTAT == 0 ) then
    KU = LDA - 2*LKL - 1

#if REALBITS==32
    call sgbsv(N,LKL,KU,NRHS,A,LDA,LPIV,B,N,LINFO)
#elif REALBITS==64
    call dgbsv(N,LKL,KU,NRHS,A,LDA,LPIV,B,N,LINFO)
#endif
  end if

  if ( .NOT.present(IPIV) ) then
    deallocate(LPIV)
  end if
end subroutine gbsv

end module lapack95
