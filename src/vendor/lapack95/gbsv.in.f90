!! ADAPTED FROM LAPACK_95 by Guy Grubbs and Michael Hirsch
module vendor_lapack95
use, intrinsic :: iso_fortran_env, only: real32, real64

implicit none (type, external)
private
public :: gbsv

interface gbsv
  procedure :: gbsv32, gbsv64
end interface gbsv

contains


subroutine gbsv32(A,B,KL,IPIV,INFO)
real(real32), dimension(:,:), intent(inout) :: A
real(real32), dimension(:), intent(inout) :: B
integer, intent(in), optional :: KL
integer, dimension(:), intent(out), optional, target :: IPIV
integer, intent(out), optional :: INFO
integer :: LINFO, ISTAT, SIPIV, LDA, N, NRHS, LKL, KU
integer, dimension(:), pointer :: LPIV
intrinsic :: size, present

external :: sgbsv

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

  @sgbsv@
endif

if(present(info)) info = linfo

if ( .NOT.present(IPIV) ) deallocate(LPIV)

end subroutine gbsv32


subroutine gbsv64(A,B,KL,IPIV,INFO)
real(real64), dimension(:,:), intent(inout) :: A
real(real64), dimension(:), intent(inout) :: B
integer, intent(in), optional :: KL
integer, dimension(:), intent(out), optional, target :: IPIV
integer, intent(out), optional :: INFO
integer :: LINFO, ISTAT, SIPIV, LDA, N, NRHS, LKL, KU
integer, dimension(:), pointer :: LPIV
intrinsic :: size, present

external :: dgbsv

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

  @dgbsv@
endif

if(present(info)) info = linfo

if ( .NOT.present(IPIV) ) deallocate(LPIV)

end subroutine gbsv64

end module vendor_lapack95
