!! ADAPTED FROM LAPACK_95 by Guy Grubbs and Michael Hirsch
module vendor_lapack95
use, intrinsic :: iso_fortran_env, only: real32, real64
implicit none

contains

subroutine gbsv(A,B,KL,IPIV,INFO)
class(*), dimension(:,:), intent(inout) :: A
class(*), dimension(:), intent(inout) :: B
integer, intent(in), optional :: KL
integer, dimension(:), intent(out), optional, target :: IPIV
integer, intent(out), optional :: INFO
integer :: LINFO, ISTAT, ISTAT1, SIPIV, LDA, N, NRHS, LKL, KU
integer, dimension(:), pointer :: LPIV
intrinsic :: size, present

LINFO = 0; ISTAT = 0;
LDA = size(A,1); N = size(A,2); NRHS = 1;

if ( present(KL) ) then
  LKL = KL
else
  LKL = (LDA-1)/3
end if

if ( present(IPIV) ) then
  SIPIV = size(IPIV)
  LPIV => IPIV
else
  SIPIV = N
  allocate(LPIV(N))
end if

if ( ISTAT == 0 ) then
  KU = LDA - 2*LKL - 1

  select type (A)
    type is (real(real32))
      call sgbsv(N,LKL,KU,NRHS,A,LDA,LPIV,B,N,LINFO)
    type is (real(real64))
      call dgbsv(N,LKL,KU,NRHS,A,LDA,LPIV,B,N,LINFO)
    class default
      error stop "unhandled kind"
  end select
endif

end subroutine gbsv

end module vendor_lapack95
