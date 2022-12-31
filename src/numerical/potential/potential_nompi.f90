!> Utility procedures not involving mpi communication (directly at least)
module potential_nompi

use phys_consts, only: wp, pi, lsp, debug, ms, qs, kB
use grid, only: gridflag, lx1,lx2,lx3,lx2all,lx3all
use meshobj, only: curvmesh
use calculus, only: grad3d2, grad3d3

implicit none (type, external)
private
public :: velocities_nompi,set_fields_test

contains
  !> This is a subroutine to compute velocities assuming that the primary state variables n,v,T have
  !    already been haloed.  
  subroutine velocities_nompi(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,flaggravdrift,flagdiamagnetic,vs2,vs3)
    !> compute steady state drifts resulting from a range of forces.  Can be used
    !   by both root and worker processes
    real(wp), dimension(:,:,:,:), intent(in) :: muP,muH,nusn
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: E2,E3
    real(wp), dimension(:,:,:), intent(in) :: vn2,vn3
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
    !! these must have ghost cells
    class(curvmesh), intent(in) :: x
    logical, intent(in) :: flaggravdrift
    logical, intent(in) :: flagdiamagnetic
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs2,vs3
    !! intent(out)
    !! these have ghost cells
    integer :: isp
    real(wp), dimension(-1:lx1+2,-1:lx2+2,-1:lx3+2) :: pressure    ! temp space for computing these
    real(wp), dimension(0:lx1+1,0:lx2+1,0:lx3+1) :: gradlp2,gradlp3

    !! For adding pressure terms we may want to just use collision freq.
    !! to avoid a bunch of slightly different mobility arrays.
    !> electric field and wind terms for ion drifts
    do isp=1,lsp
      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2(1:lx1,1:lx2,1:lx3)-muH(:,:,:,isp)*E3(1:lx1,1:lx2,1:lx3)+ &
                        (muP(:,:,:,isp)*vn2-muH(:,:,:,isp)*vn3)*(ms(isp)*nusn(:,:,:,isp)/qs(isp))
      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2(1:lx1,1:lx2,1:lx3)+muP(:,:,:,isp)*E3(1:lx1,1:lx2,1:lx3)+ &
                        (muH(:,:,:,isp)*vn2+muP(:,:,:,isp)*vn3)*ms(isp)*nusn(:,:,:,isp)/qs(isp)
    end do

    !> Pressure/diamagnetic terms (if required)
    if (flagdiamagnetic) then
      do isp=1,lsp
        !> this behaves better when we take the gradient of log pressure
        pressure(-1:lx1+2,-1:lx2+2,-1:lx3+2)=log(ns(-1:lx1+2,-1:lx2+2,-1:lx3+2,isp)*kB*Ts(-1:lx1+2,-1:lx2+2,-1:lx3+2,isp))
        gradlp2=grad3D2(pressure(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
        gradlp3=grad3D3(pressure(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
        vs2(1:lx1,1:lx2,1:lx3,isp)=vs2(1:lx1,1:lx2,1:lx3,isp) &
                   -muP(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradlp2(1:lx1,1:lx2,1:lx3) &
                   +muH(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradlp3(1:lx1,1:lx2,1:lx3)
        vs3(1:lx1,1:lx2,1:lx3,isp)=vs3(1:lx1,1:lx2,1:lx3,isp) &
                   -muH(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradlp2(1:lx1,1:lx2,1:lx3) &
                   -muP(1:lx1,1:lx2,1:lx3,isp)*kB*Ts(1:lx1,1:lx2,1:lx3,isp)/qs(isp)*gradlp3(1:lx1,1:lx2,1:lx3)
      end do
    end if

    !> Gravitational drift terms (if required)
    if (flaggravdrift) then
      do isp=1,lsp
        vs2(1:lx1,1:lx2,1:lx3,isp)=vs2(1:lx1,1:lx2,1:lx3,isp)+ms(isp)/qs(isp)*(muP(:,:,:,isp)*x%g2-muH(:,:,:,isp)*x%g3)    !FIXME: +muH looks suspicious, I'm changing to (-)
        vs3(1:lx1,1:lx2,1:lx3,isp)=vs3(1:lx1,1:lx2,1:lx3,isp)+ms(isp)/qs(isp)*(muH(:,:,:,isp)*x%g2+muP(:,:,:,isp)*x%g3)
      end do
    end if


    !! If it were appropriate this is how polarzations drifts could be computed.  However the particular quasistatic
    !   model that we use explicitly omits this from the drift calculation which is then used in convective term in
    !   polarization current.  Physically it accounts for charge accumulation from polarization currents but not for
    !   the *direct* effect of the polarization term on drift.
    !    do isp=1,lsp
           !! To leading order the ion drifts do not include the polarization parts,
           !! otherwise it may mess up polarization convective term in the electrodynamics solver...
    !      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE2Dt
    !      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE3Dt
    !    end do
  end subroutine velocities_nompi


  !> Set the electric fields to some fixed value ***for purposes of testing***; this shouldn't be used for any other purpose
  subroutine set_fields_test(x,E1,E2,E3)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: E1,E2,E3
    integer ix1,ix2,ix3
    real(wp) :: r0=50e3,sigr=10e3
    real(wp) :: Er,r
    real(wp) :: Eamp=50e-3
    real(wp) :: phiang

    E1=0._wp
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          r=sqrt(x%x2(ix2)**2+x%x3(ix3)**2)
          phiang=atan2(x%x3(ix3),x%x2(ix2))
          Er=Eamp*exp(-(r-r0)**2/2/sigr**2)
          if (r<10e3 .or. r>90e3) Er=0._wp

          E2(ix1,ix2,ix3)=Er*cos(phiang)
          E3(ix1,ix2,ix3)=Er*sin(phiang)
        end do
      end do
    end do
  end subroutine set_fields_test
end module potential_nompi
