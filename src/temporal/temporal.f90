module temporal

use phys_consts, only: kB,mu0,ms,lsp,pi, wp, debug
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

implicit none (type, external)

private
public :: cflcalc

contains
  !> Compute the max cfl number of the entirety of the worker grid
  subroutine cflcalc(Ts,vs1,vs2,vs3,dl1i,dl2i,dl3i,dt,maxcfl)
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: Ts,vs1,vs2,vs3
    real(wp), dimension(:,:,:), intent(in) :: dl1i
    real(wp), dimension(:,:,:), intent(in) :: dl2i
    real(wp), dimension(:,:,:), intent(in) :: dl3i
    real(wp), intent(in) :: dt
    real(wp), intent(out) :: maxcfl
    real(wp) :: vsnd
    integer :: lx1,lx2,lx3,ix1,ix2,ix3,isp
    real(wp) :: cfltmp

    lx1=size(Ts,1)-4
    lx2=size(Ts,2)-4
    lx3=size(Ts,3)-4

    if (.not.ieee_is_finite(dt)) error stop "cflcalc: nonfinite step"
    if (dt<=0) error stop "cflcalc: nonpositive step"
    if (min(lx1,lx2,lx3)<1 .or. size(Ts,4)/=lsp) error stop "cflcalc: invalid state shape"
    if (any(shape(Ts)/=shape(vs1)) .or. any(shape(Ts)/=shape(vs2)) .or. &
        any(shape(Ts)/=shape(vs3))) error stop "cflcalc: velocity shape mismatch"
    if (any(shape(dl1i)/=[lx1,lx2,lx3]) .or. any(shape(dl2i)/=[lx1,lx2,lx3]) .or. &
        any(shape(dl3i)/=[lx1,lx2,lx3])) error stop "cflcalc: metric shape mismatch"
    if (.not.all(ieee_is_finite(Ts(1:lx1,1:lx2,1:lx3,:))) .or. &
        .not.all(ieee_is_finite(vs1(1:lx1,1:lx2,1:lx3,:))) .or. &
        .not.all(ieee_is_finite(vs2(1:lx1,1:lx2,1:lx3,:))) .or. &
        .not.all(ieee_is_finite(vs3(1:lx1,1:lx2,1:lx3,:)))) error stop "cflcalc: nonfinite state"
    if (minval(Ts(1:lx1,1:lx2,1:lx3,:))<0) error stop "cflcalc: negative temperature"
    if (.not.all(ieee_is_finite(dl1i)) .or. .not.all(ieee_is_finite(dl2i)) .or. &
        .not.all(ieee_is_finite(dl3i))) error stop "cflcalc: nonfinite metric"
    if (min(minval(dl1i),minval(dl2i),minval(dl3i))<=0) error stop "cflcalc: nonpositive metric"

    !EVALUATE TIME STEP AGAINST LOCAL SOUND SPEED AND ADVECTION
    maxcfl=0._wp
    do isp=1,lsp
      do ix3=1,lx3
        do ix2=1,lx2
          do ix1=1,lx1
            if (isp<lsp) then
              vsnd=sqrt(kB*Ts(ix1,ix2,ix3,isp)/ms(isp)+5._wp/3._wp*kB*Ts(ix1,ix2,ix3,lsp)/ms(isp))
            else
              vsnd=0._wp
            end if

            cfltmp=dt*(vsnd+abs(vs1(ix1,ix2,ix3,isp)))/dl1i(ix1,ix2,ix3)
            if (cfltmp>maxcfl) maxcfl=cfltmp
            cfltmp=dt*abs(vs2(ix1,ix2,ix3,isp))/dl2i(ix1,ix2,ix3)
            if (cfltmp>maxcfl) maxcfl=cfltmp
            cfltmp=dt*abs(vs3(ix1,ix2,ix3,isp))/dl3i(ix1,ix2,ix3)
            if (cfltmp>maxcfl) maxcfl=cfltmp
          end do
        end do
      end do
    end do
  end subroutine cflcalc
end module temporal
