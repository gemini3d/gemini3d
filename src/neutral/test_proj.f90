!> tests the projections for a given grid structure
program test_proj

use pathlib, only : mkdir
use phys_consts, only: wp
use meshobj_dipole, only : dipolemesh
use neutral, only : store_geo2native_projections

implicit none (type, external)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, parameter :: lq=384+4,lp=96+4,lphi=64+4
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(lphi) :: phi
! from tohoku20113D_lowres_3Dneu
real(wp), dimension(2), parameter :: qlims=[-0.5340405,0.5340405]
real(wp), dimension(2), parameter :: plims=[1.2509838,1.4372374]
real(wp), dimension(2), parameter :: philims=[3.6126509,3.7240195]
integer :: iq,ip,iphi, i
real(wp) :: minchkvar,maxchkvar
real(wp), dimension(:,:,:), allocatable :: proj
character(:), allocatable :: path
character(1000) :: argv
type(dipolemesh) :: x
real(wp), dimension(lq,lp,lphi,3) :: ealt,eglat,eglon
real(wp), dimension(3,3,lq,lp,lphi) :: rotmats
real(wp), dimension(3,3) :: matnow,eyetest,eye
real(wp), parameter :: errthresh=1e-6
real(wp) :: eyeerr
logical :: debug=.false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! create a dipole grid
q=[(qlims(1) + (qlims(2)-qlims(1))/(lq-1)*(iq-1),iq=1,lq)]
p=[(plims(1) + (plims(2)-plims(1))/(lp-1)*(ip-1),ip=1,lp)]
phi=[(philims(1) + (philims(2)-philims(1))/(lphi-1)*(iphi-1),iphi=1,lphi)]
call x%set_coords(q,p,phi,p,phi)
call x%init()
call x%make()

!! compute the geographic projections
call x%calc_unitvec_geo(ealt,eglon,eglat)
call store_geo2native_projections(x,ealt,eglon,eglat,rotmat=rotmats)

!! verify that the transformation is approximately unitary at all non-ghost grid points
print*, 'Begin testing transformation...'
eye=reshape([1,0,0,0,1,0,0,0,1],[3, 3])
do iphi=1,lphi-4
  do ip=1,lp-4
    do iq=1,lq-4
      matnow=rotmats(1:3,1:3,iq,ip,iphi)
      eyetest=matmul(matnow,transpose(matnow))     ! rotations are unitary so R*R.transpose() should be identity...
      eyeerr=sum(abs(eyetest-eye))                 ! aggregate error evaluation
      if (debug) call printmats()
      if (eyeerr > errthresh) then
        print*, 'Excessive deviation from unitary:  ',eyeerr
        print*, '  at grid indices:  ',iq,ip,iphi
        print*, '-------------------------------------------------------'
        call printmats()
        error stop
      end if
    end do
  end do
end do
print*, 'Transformation tests complete successfully...'

contains
  subroutine printmats()
    print*, eye(1,1:3)
    print*, eye(2,1:3)
    print*, eye(3,1:3)
    print*, '-------------------------------------------------------'
    print*, eyetest(1,1:3)
    print*, eyetest(2,1:3)
    print*, eyetest(3,1:3)
    print*, '-------------------------------------------------------'
  end subroutine printmats
end program test_proj

