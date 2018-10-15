!------ ADAPTED FROM LAPACK_95 by Guy Grubbs --------
module f95_lapack
use phys_consts, only: wp
implicit none
contains

subroutine la_gbsv(A,B,KL,IPIV,INFO)
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
    call dgbsv(N,LKL,KU,NRHS,A,LDA,LPIV,B,N,LINFO)
  end if

  if ( .NOT.present(IPIV) ) then
    deallocate(LPIV)
  end if
end subroutine la_gbsv

end module f95_lapack
