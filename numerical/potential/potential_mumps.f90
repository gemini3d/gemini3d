module potential_mumps


!NOTES:
! - IF ONE REALLY WANTED TO CLEAN THIS UP IT MIGHT BE MORE EFFICIENT TO USE HARWELL-BOEING FORMAT
!   FOR MATRICES...

!SOME SUPERFLUOUS ARGUMENTS THAT ARE LEFT IN TO MAINTAIN UNIFORMITY ACROSS CALLS
!/home/zettergm/zettergmdata/GEMINI/numerical/potential/potential_mumps.f90:1291:0:
!warning: unused parameter ‘vminx1’ [-Wunused-parameter]
!   function
!elliptic2D_nonint_curv(srcterm,sig0,sigP,Vminx1,Vmaxx1,Vminx3,Vmaxx3,x,flagdirich,perflag,it)
! ^
!/home/zettergm/zettergmdata/GEMINI/numerical/potential/potential_mumps.f90:764:0:
!warning: unused parameter ‘vminx3’ [-Wunused-parameter]
!   function
!elliptic2D_pol_conv_curv_periodic2(srcterm,SigP,SigH,Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x,Phi0,perflag,it)
! ^
!/home/zettergm/zettergmdata/GEMINI/numerical/potential/potential_mumps.f90:764:0:
!warning: unused parameter ‘vmaxx3’ [-Wunused-parameter]
!/home/zettergm/zettergmdata/GEMINI/numerical/potential/potential_mumps.f90:22:0:
!warning: unused parameter ‘vminx1’ [-Wunused-parameter]
!   function
!elliptic3D_curv(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3,
!&
! ^

use, intrinsic:: iso_fortran_env, only: stderr=>error_unit, stdout=>output_unit

use mpi, only: mpi_comm_world
use phys_consts, only: wp, debug
use calculus, only: grad3D1, grad3D2, grad3D3, grad2D1_curv_alt, grad2D3, grad2D3_curv_periodic
use grid, only: gridflag
use mesh, only: curvmesh
use mpimod, only: myid
use interpolation, only: interp1
use PDEelliptic, only: elliptic3D_cart,elliptic2D_polarization,elliptic2D_polarization_periodic,elliptic2D_cart

implicit none

private

integer, dimension(:), pointer, protected, save :: mumps_perm   !cached permutation, unclear whether save is necessary...

public :: potential3D_fieldresolved_decimate, potential2D_polarization, potential2D_polarization_periodic, potential2D_fieldresolved

contains


function potential3D_fieldresolved_decimate(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                  x,flagdirich,perflag,it)

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
!-------ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
!-------LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
!-------GRID POINTS.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP,sigH
real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
type(curvmesh), intent(in) :: x
integer, intent(in) :: flagdirich
logical, intent(in) :: perflag
integer, intent(in) :: it

real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigP2,gradsigP3
real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigH2,gradsigH3
real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsig01
real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: Ac,Bc,Cc,Dc,Ec,Fc

integer :: lx1,lx2,lx3,ix1,ix2,ix3
integer, parameter :: ldec=11
real(wp), dimension(:), allocatable :: x1dec
real(wp), dimension(:), allocatable :: dx1dec
real(wp), dimension(:), allocatable :: x1idec
real(wp), dimension(:), allocatable :: dx1idec
real(wp), dimension(:,:,:), allocatable :: Acdec,Bcdec,Ccdec,Dcdec,Ecdec,Fcdec,srctermdec
real(wp), dimension(:,:), allocatable :: Vminx2dec,Vmaxx2dec
real(wp), dimension(:,:), allocatable :: Vminx3dec, Vmaxx3dec
real(wp), dimension(:,:,:), allocatable :: Phidec

real(wp), dimension(1:size(Vminx1,1),1:size(Vminx1,2)) :: Vminx1pot,Vmaxx1pot

real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: potential3D_fieldresolved_decimate


!SYSTEM SIZES
lx1=x%lx1    !These will be full grid sizes if called from root (only acceptable thing)
lx2=x%lx2all
lx3=x%lx3all


!COMPUTE AUXILIARY COEFFICIENTS TO PASS TO CART SOLVER
if (debug) print *, 'Prepping coefficients for elliptic equation...'
gradsig01=grad3D1(sig0,x,1,lx1,1,lx2,1,lx3)
gradsigP2=grad3D2(sigP,x,1,lx1,1,lx2,1,lx3)
gradsigP3=grad3D3(sigP,x,1,lx1,1,lx2,1,lx3)
gradsigH2=grad3D2(sigH,x,1,lx1,1,lx2,1,lx3)
gradsigH3=grad3D3(sigH,x,1,lx1,1,lx2,1,lx3)

Ac=sigP
Bc=sigP
Cc=sig0
Dc=gradsigP2+gradsigH3
Ec=gradsigP3-gradsigH2
Fc=gradsig01


!DEFINE A DECIMATED MESH (THIS IS HARDCODED FOR NOW)
if (debug) print*, 'Decimating parallel grid...'
allocate(x1dec(-1:ldec+2),dx1dec(0:ldec+2),x1idec(1:ldec+1),dx1idec(1:ldec))
x1dec(-1:lx1+2)=[x%x1(-1),x%x1(0),x%x1(1),81.8e3_wp,84.2e3_wp,87.5e3_wp,93.3e3_wp,106.0e3_wp,124.0e3_wp, &
            144.6e3_wp,206.7e3_wp,882.2e3_wp,x%x1(lx1),x%x1(lx1+1),x%x1(lx1+2)]
!x1dec(-1:lx1+2)=[x%x1(-1),x%x1(0),x%x1(1),81.8e3_wp,84.2e3_wp,87.5e3_wp,93.3e3_wp,106.0e3_wp,124.0e3_wp, &
!            144.6e3_wp,175e3_wp,206.7e3_wp,250e3_wp,400e3_wp,600e3_wp,882.2e3_wp,x%x1(lx1),x%x1(lx1+1),x%x1(lx1+2)]

dx1dec(0:ldec+2)=x1dec(0:ldec+2)-x1dec(-1:ldec+1)
x1idec(1:ldec+1)=0.5_wp*(x1dec(0:ldec)+x1dec(1:ldec+1))
dx1idec(1:ldec)=x1idec(2:ldec+1)-x1idec(1:ldec)


!INTERPOLATE COEFFICIENTS AND SOURCE TERM ONTO DECIMATED GRID
if (debug) print*, 'Interpolating coefficients...'
allocate(Acdec(1:ldec,1:lx2,1:lx3),Bcdec(1:ldec,1:lx2,1:lx3),Ccdec(1:ldec,1:lx2,1:lx3), &
         Dcdec(1:ldec,1:lx2,1:lx3),Ecdec(1:ldec,1:lx2,1:lx3),Fcdec(1:ldec,1:lx2,1:lx3), &
         srctermdec(1:ldec,1:lx2,1:lx3))
do ix2=1,lx2
  do ix3=1,lx3
    Acdec(:,ix2,ix3)=interp1(x%x1(1:lx1),Ac(:,ix2,ix3),x1dec(1:ldec))
    Bcdec(:,ix2,ix3)=interp1(x%x1(1:lx1),Bc(:,ix2,ix3),x1dec(1:ldec))
    Ccdec(:,ix2,ix3)=interp1(x%x1(1:lx1),Cc(:,ix2,ix3),x1dec(1:ldec))
    Dcdec(:,ix2,ix3)=interp1(x%x1(1:lx1),Dc(:,ix2,ix3),x1dec(1:ldec))
    Ecdec(:,ix2,ix3)=interp1(x%x1(1:lx1),Ec(:,ix2,ix3),x1dec(1:ldec))
    Fcdec(:,ix2,ix3)=interp1(x%x1(1:lx1),Fc(:,ix2,ix3),x1dec(1:ldec))
    srctermdec(:,ix2,ix3)=interp1(x%x1(1:lx1),srcterm(:,ix2,ix3),x1dec(1:ldec))
  end do
end do


!INTERPOLATE BOUNDARY CONDITIONS ONTO DECIMATED GRID
allocate(Vminx2dec(1:ldec,1:lx3),Vmaxx2dec(1:ldec,1:lx3))
do ix3=1,lx3
  Vminx2dec(:,ix3)=interp1(x%x1(1:lx1),Vminx2(:,ix3),x1dec(1:ldec))
  Vmaxx2dec(:,ix3)=interp1(x%x1(1:lx1),Vmaxx2(:,ix3),x1dec(1:ldec))
end do
allocate(Vminx3dec(1:ldec,1:lx2),Vmaxx3dec(1:ldec,1:lx2))
do ix2=1,lx2
  Vminx3dec(:,ix2)=interp1(x%x1(1:lx1),Vminx3(:,ix2),x1dec(1:ldec))
  Vmaxx3dec(:,ix2)=interp1(x%x1(1:lx1),Vmaxx3(:,ix2),x1dec(1:ldec))
end do


!FOR WHATEVER REASON THE EDGE VALUES GET MESSED UP BY INTERP1
Acdec(1,:,:)=Ac(1,:,:)
Acdec(ldec,:,:)=Ac(lx1,:,:)
Bcdec(1,:,:)=Bc(1,:,:)
Bcdec(ldec,:,:)=Bc(lx1,:,:)
Ccdec(1,:,:)=Cc(1,:,:)
Ccdec(ldec,:,:)=Cc(lx1,:,:)
Dcdec(1,:,:)=Dc(1,:,:)
Dcdec(ldec,:,:)=Dc(lx1,:,:)
Ecdec(1,:,:)=Ec(1,:,:)
Ecdec(ldec,:,:)=Ec(lx1,:,:)
Fcdec(1,:,:)=Fc(1,:,:)
Fcdec(ldec,:,:)=Fc(lx1,:,:)
Vminx2dec(1,:)=Vminx2(1,:)
Vminx2dec(ldec,:)=Vminx2(lx1,:)
Vmaxx2dec(1,:)=Vmaxx2(1,:)
Vmaxx2dec(ldec,:)=Vmaxx2(lx1,:)
Vminx3dec(1,:)=Vminx3(1,:)
Vminx3dec(ldec,:)=Vminx3(lx1,:)
Vmaxx3dec(1,:)=Vmaxx3(1,:)
Vmaxx3dec(ldec,:)=Vmaxx3(lx1,:)
srctermdec(1,:,:)=srcterm(1,:,:)
srctermdec(ldec,:,:)=srcterm(lx1,:,:)

!
!print*, minval(Acdec),maxval(Acdec)
!print*, minval(Ac),maxval(Ac)
!print*, minval(Bcdec),maxval(Bcdec)
!print*, minval(Bc),maxval(Bc)
!print*, minval(Ccdec),maxval(Ccdec)
!print*, minval(Cc),maxval(Cc)
!print*, minval(Dcdec),maxval(Dcdec)
!print*, minval(Dc),maxval(Dc)
!print*, minval(Ecdec),maxval(Ecdec)
!print*, minval(Ec),maxval(Ec)
!print*, minval(Fcdec),maxval(Fcdec)
!print*, minval(Fc),maxval(Fc)
!print*, minval(srctermdec),maxval(srctermdec)
!print*, minval(srcterm),maxval(srcterm)
!print*, minval(Vminx2dec),maxval(Vminx2dec)
!print*, minval(Vmaxx2dec),maxval(Vmaxx2dec)
!print*, minval(Vminx3dec),maxval(Vminx3dec)
!print*, minval(Vmaxx3dec),maxval(Vmaxx3dec)
!print*, minval(dx1dec),maxval(dx1dec)
!print*, minval(dx1idec),maxval(dx1idec)
!print*, x1dec(-1:ldec+2)
!print*, dx1dec(0:ldec+2)
!

!ADJUST THE BOUNDARY CONDITION TO POTENTIAL DERIVATIVE INSTEAD OF CURRENT DENSITY
Vminx1pot=-1*Vminx1/sig0(1,:,:)
Vmaxx1pot=-1*Vmaxx1/sig0(lx1,:,:)


!CALL CARTESIAN SOLVER ON THE DECIMATED GRID
if (debug) print*, 'Calling solve on decimated grid...'
allocate(Phidec(1:ldec,1:lx2,1:lx3))
Phidec=elliptic3D_cart(srctermdec,Acdec,Bcdec,Ccdec,Dcdec,Ecdec,Fcdec,Vminx1pot,Vmaxx1pot, &
                Vminx2dec,Vmaxx2dec,Vminx3dec,Vmaxx3dec, &
                dx1dec,dx1idec,x%dx2all,x%dx2iall,x%dx3all,x%dx3iall,flagdirich,perflag,it)


!INTERPOLATE BACK UP TO MAIN GRID
if (debug) print*, 'Upsampling potential...'
do ix2=1,lx2
  do ix3=1,lx3
    potential3D_fieldresolved_decimate(:,ix2,ix3)=interp1(x1dec(1:ldec),Phidec(:,ix2,ix3),x%x1(1:lx1))
  end do
end do


!AGAIN NEED TO FIX THE EDGES...
potential3D_fieldresolved_decimate(1,:,:)=Phidec(1,:,:)
potential3D_fieldresolved_decimate(lx1,:,:)=Phidec(ldec,:,:)


! open(newunit=u, form='unformatted', access='stream',file='Phidec.raw8',status='replace', action='write')
! write(u) potential3D_fieldresolved_decimate,Phidec,Ac,Acdec,Bc,Bcdec,Cc,Ccdec,Dc,Dcdec,Ec,Ecdec,Fc,Fcdec,srcterm,srctermdec
! write(u) Vminx1pot,Vmaxx1pot,Vminx2dec,Vmaxx2dec,Vminx3dec,Vmaxx3dec
! close(u)


!CLEAN UP THE ALLOCATED ARRAYS
deallocate(srctermdec,Acdec,Bcdec,Ccdec,Dcdec,Ecdec,Fcdec,dx1dec,x1dec,x1idec,dx1idec,Vminx2dec,Vmaxx2dec,Vminx3dec,Vmaxx3dec)
deallocate(Phidec)

end function potential3D_fieldresolved_decimate


function potential2D_polarization(srcterm,SigP2,SigP3,SigH,Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x,Phi0,perflag,it)

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS
!-------INCLUDES FULL OF POLARIZATION CURRENT, INCLUDING CONVECTIVE
!-------TERMS.  VELOCITIES SHOULD BE TRIMMED (WITHOUT GHOST CELLS).
!-------THIS VERSION OF THE *INTEGRATED* POTENTIAL SOLVER OBVIATES
!-------ALL OTHERS SINCE A PURELY ELECTRSTATIC FORM CAN BE RECOVERED
!-------BY ZEROING OUT THE INERTIAL CAPACITANCE.
!-------
!-------THIS FORM IS INTENDED TO  WORK WITH CURVILINEAR MESHES.
!-------NOTE THAT THE FULL GRID VARIABLES (X%DX3ALL, ETC.) MUST
!-------BE USED HERE!!!
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: srcterm,SigP2,SigP3,SigH,Cm,v2,v3    !ZZZ - THESE WILL NEED TO BE MODIFIED CONDUCTIVITIES, AND WE'LL NEED THREE OF THEM
real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
real(wp), dimension(:), intent(in) :: Vminx3,Vmaxx3
real(wp), intent(in) :: dt
type(curvmesh), intent(in) :: x
real(wp), dimension(:,:), intent(in) :: Phi0
logical, intent(in) :: perflag
integer, intent(in) :: it

real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: gradSigH2,gradSigH3
integer :: u
integer :: lx2,lx3

real(wp), dimension(size(SigP2,1),size(SigP2,2)) :: potential2D_polarization


lx2=x%lx2all    !use full grid sizes
lx3=x%lx3all

!gradSigH2=grad2D1(SigH,x,1,lx2)   !x2 is now 1st index and x3 is second...  This one appears to be the problem.  This issue here is that grad2D1 automatically uses x%dx1 as the differential element...
gradSigH2=grad2D1_curv_alt(SigH,x,1,lx2)
!! note the alt since we need to use dx2 as differential...  Tricky bug/feature
gradSigH3=grad2D3(SigH,x,1,lx3)
!! awkward way of handling this special case derivative which uses x3 as the differential to operate on a 2D array.

potential2D_polarization=elliptic2D_polarization(srcterm,SigP2,SigP3,SigH,gradSigH2,gradSigH3,Cm,v2,v3, &
                                  Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x%dx1,x%dx1i,x%dx2all,x%dx2iall, &
                                  x%dx3all,x%dx3iall,Phi0,perflag,it)

end function potential2D_polarization


function potential2D_polarization_periodic(srcterm,SigP,SigH,Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x,Phi0,perflag,it)

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS
!-------INCLUDES FULL OF POLARIZATION CURRENT, INCLUDING CONVECTIVE
!-------TERMS.  VELOCITIES SHOULD BE TRIMMED (WITHOUT GHOST CELLS).
!-------THIS VERSION OF THE *INTEGRATED* POTENTIAL SOLVER OBVIATES
!-------ALL OTHERS SINCE A PURELY ELECTRSTATIC FORM CAN BE RECOVERED
!-------BY ZEROING OUT THE INERTIAL CAPACITANCE.
!-------
!-------THIS FORM IS INTENDED TO  WORK WITH CARTESIAN MESHES ONLY.
!-------NOTE THAT THE FULL GRID VARIABLES (X%DX3ALL, ETC.) MUST
!-------BE USED HERE!!!
!-------
!-------THIS FUNCTION WORKS ON A PERIODIC MESH BY USING A CIRCULANT MATRIX
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: srcterm,SigP,SigH,Cm,v2,v3
real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
real(wp), dimension(:), intent(in) :: Vminx3,Vmaxx3
real(wp), intent(in) :: dt
type(curvmesh), intent(in) :: x
real(wp), dimension(:,:), intent(in) :: Phi0
logical, intent(in) :: perflag
integer, intent(in) :: it

real(wp), dimension(1:size(SigP,1),1:size(SigP,2)+1) :: gradSigH2,gradSigH3
real(wp), dimension(size(SigP,1),size(SigP,2)) :: potential2D_polarization_periodic
integer :: lx2,lx3


lx2=x%lx2all    !use full grid sizes
lx3=x%lx3all

!ZZZ - THESE NEED TO BE CHANGED INTO CIRCULAR/PERIODIC DERIVATIVES FOR THE X3 DIRECTION
gradSigH2=grad2D1_curv_alt(SigH,x,1,lx2)   !note the alt since we need to use dx2 as differential...  Tricky bug/feature
gradSigH3=grad2D3_curv_periodic(SigH,x,1,lx3)    !circular difference


potential2D_polarization_periodic=elliptic2D_polarization_periodic(srcterm,SigP,SigH,gradSigH2,gradSigH3, &
                          Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x%dx1,x%dx1i,x%dx2all,x%dx2iall,x%dx3all, &
                          x%dx3iall,Phi0,perflag,it)

end function potential2D_polarization_periodic


function potential2D_fieldresolved(srcterm,sig0,sigP,Vminx1,Vmaxx1,Vminx3,Vmaxx3,x,flagdirich,perflag,it)

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS
!-------ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
!-------LINE AND THAT IT VARIES IN X1 AND X3 (X2 IS NOMINALL JUST
!-------ONE ELEMENT.  LEFT AND RIGHT BOUNDARIES (IN X3) ARE ASSUMED
!-------TO USE DIRICHLET BOUNARY CONDITIONS, WHILE THE (ALTITUDE)
!-------TOP CAN BE NEUMANN OR DIRICHLET.  BOTTOM (ALTITUDE)
!-------IS ALWAYS ASSUMED TO BE DIRICHLET.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP   !arrays passed in will still have full rank 3
real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
type(curvmesh), intent(in) :: x
integer, intent(in) :: flagdirich
logical, intent(in) :: perflag
integer, intent(in) :: it

real(wp), dimension(size(Vmaxx1,1),size(Vmaxx1,2)) :: Vmaxx1alt      !in case we need to do some transformations to convert current into potential
real(wp), dimension(size(sig0,1),1,size(sig0,3)) :: potential2D_fieldresolved
integer :: lx1

lx1=size(sig0,1)
if (flagdirich==0) then      !convert current into potential normal derivative
  Vmaxx1alt=-1*Vmaxx1*x%h1(lx1,:,:)/sig0(lx1,:,:)
else                         !Dirichlet boundary conditions, don't change
  Vmaxx1alt=Vmaxx1
end if

potential2D_fieldresolved=elliptic2D_cart(srcterm,sig0,sigP,Vminx1,Vmaxx1alt,Vminx3,Vmaxx3, &
                     x%dx1,x%dx1i,x%dx3all,x%dx3iall,flagdirich,perflag,gridflag,it)

end function potential2D_fieldresolved

end module potential_mumps
