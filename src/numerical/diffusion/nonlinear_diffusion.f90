module nonlinear_diffusion
! Proposed nonlinear time centering. The predictor is virtual: only the final
! TRBDF2 increments describe applied changes and enter the numerical ledger.
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp
use meshobj, only: curvmesh
use collisions, only: thermal_conduct
use diffusion, only: diffusion_prep, backEuler3D, TRBDF23D
implicit none(type,external)
private
public :: midpoint_TRBDF23D
contains
function midpoint_TRBDF23D(isp,x,ns,T,nn,J1,Tn,Teinf,dt,increments,lambda_applied) result(Tnew)
  integer, intent(in) :: isp
  class(curvmesh), intent(in) :: x
  real(wp), intent(in) :: ns(-1:,-1:,-1:),T(-1:,-1:,-1:),nn(:,:,:,:),J1(-1:,-1:,-1:),Tn(:,:,:)
  real(wp), intent(in) :: Teinf,dt
  real(wp), intent(out), optional :: increments(size(T,1)-4,size(T,2)-4,size(T,3)-4,5)
  real(wp), intent(out), optional :: lambda_applied(size(T,1)-4,size(T,2)-4,size(T,3)-4)
  real(wp) :: Tnew(-1:size(T,1)-2,-1:size(T,2)-2,-1:size(T,3)-2)
  real(wp), allocatable :: initial(:,:,:),predictor(:,:,:),lambda(:,:,:),beta(:,:,:)
  real(wp), allocatable :: A(:,:,:),B(:,:,:),C(:,:,:),D(:,:,:),E(:,:,:)
  integer :: n1,n2,n3
  n1=size(T,1)-4;n2=size(T,2)-4;n3=size(T,3)-4
  if(.not.ieee_is_finite(dt).or.dt<=0) error stop 'midpoint diffusion requires positive finite dt'
  allocate(initial(-1:n1+2,-1:n2+2,-1:n3+2),predictor(-1:n1+2,-1:n2+2,-1:n3+2))
  allocate(lambda(n1,n2,n3),beta(n1,n2,n3),A(n1,n2,n3),B(n1,n2,n3),C(n1,n2,n3),D(n1,n2,n3),E(n1,n2,n3))
  initial=T
  call thermal_conduct(isp,initial,ns,nn,J1,lambda,beta)
  call diffusion_prep(isp,x,lambda,beta,ns,initial,A,B,C,D,E,Tn,Teinf)
  predictor=backEuler3D(initial,A,B,C,D,E,dt/2,x)
  ! Do not quietly introduce a new temperature floor in the virtual stage.
  if(any(.not.ieee_is_finite(predictor(1:n1,1:n2,1:n3))).or. &
     any(predictor(1:n1,1:n2,1:n3)<=0)) error stop 'nonpositive or nonfinite midpoint temperature'
  call thermal_conduct(isp,predictor,ns,nn,J1,lambda,beta)
  if(any(.not.ieee_is_finite(lambda)).or.any(lambda<=0)) error stop 'invalid midpoint conductivity'
  ! Rebuild geometry and conductivity-dependent Neumann ghost differences.
  ! Start the applied step from the original temperature, not the predictor.
  call diffusion_prep(isp,x,lambda,beta,ns,initial,A,B,C,D,E,Tn,Teinf)
  if(present(increments)) then
    Tnew=TRBDF23D(initial,A,B,C,D,E,dt,x,increments)
  else
    Tnew=TRBDF23D(initial,A,B,C,D,E,dt,x)
  endif
  if(present(lambda_applied)) lambda_applied=lambda
end function midpoint_TRBDF23D
end module nonlinear_diffusion
