module temporal

use phys_consts, only: kB,mu0,ms,lsp,pi, wp, debug

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
