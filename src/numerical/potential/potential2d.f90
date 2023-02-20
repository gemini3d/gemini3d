submodule (potential_mumps) potential2d

use grid, only: gridflag
use calculus, only : grad2D1_curv_alt, grad2D3, grad2D3_curv_periodic
use PDEelliptic, only: elliptic2D_polarization,elliptic2D_polarization_periodic,elliptic2D_cart

implicit none (type, external)

contains


module procedure potential2D_polarization
!! SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS
!! INCLUDES FULL OF POLARIZATION CURRENT, INCLUDING CONVECTIVE
!! TERMS.  VELOCITIES SHOULD BE TRIMMED (WITHOUT GHOST CELLS).
!! THIS VERSION OF THE *INTEGRATED* POTENTIAL SOLVER OBVIATES
!! ALL OTHERS SINCE A PURELY ELECTRSTATIC FORM CAN BE RECOVERED
!! BY ZEROING OUT THE INERTIAL CAPACITANCE.
!!
!! THIS FORM IS INTENDED TO  WORK WITH CURVILINEAR MESHES.
!! NOTE THAT THE FULL GRID VARIABLES (X%DX3ALL, ETC.) MUST
!! BE USED HERE!!!

real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: gradSigH2,gradSigH3
integer :: lx2,lx3


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

end procedure potential2D_polarization


module procedure potential2D_polarization_periodic

!! SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS
!! INCLUDES FULL OF POLARIZATION CURRENT, INCLUDING CONVECTIVE
!! TERMS.  VELOCITIES SHOULD BE TRIMMED (WITHOUT GHOST CELLS).
!! THIS VERSION OF THE *INTEGRATED* POTENTIAL SOLVER OBVIATES
!! ALL OTHERS SINCE A PURELY ELECTRSTATIC FORM CAN BE RECOVERED
!! BY ZEROING OUT THE INERTIAL CAPACITANCE.
!!
!! THIS FORM IS INTENDED TO  WORK WITH CARTESIAN MESHES ONLY.
!! NOTE THAT THE FULL GRID VARIABLES (X%DX3ALL, ETC.) MUST
!! BE USED HERE!!!
!!
!! THIS FUNCTION WORKS ON A PERIODIC MESH BY USING A CIRCULANT MATRIX

real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: gradSigH2,gradSigH3

integer :: lx2,lx3


lx2=x%lx2all    !use full grid sizes
lx3=x%lx3all

!ZZZ - THESE NEED TO BE CHANGED INTO CIRCULAR/PERIODIC DERIVATIVES FOR THE X3 DIRECTION
gradSigH2=grad2D1_curv_alt(SigH,x,1,lx2)   !note the alt since we need to use dx2 as differential...  Tricky bug/feature
gradSigH3=grad2D3_curv_periodic(SigH,x,1,lx3)    !circular difference


potential2D_polarization_periodic=elliptic2D_polarization_periodic(srcterm,SigP,SigH,gradSigH2,gradSigH3, &
                          Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x%dx1,x%dx1i,x%dx2all,x%dx2iall,x%dx3all, &
                          x%dx3iall,Phi0,perflag,it)

end procedure potential2D_polarization_periodic


module procedure potential2D_fieldresolved
!! SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS
!! ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
!! LINE AND THAT IT VARIES IN X1 AND X3 (X2 IS NOMINALL JUST
!! ONE ELEMENT.  LEFT AND RIGHT BOUNDARIES (IN X3) ARE ASSUMED
!! TO USE DIRICHLET BOUNDARY CONDITIONS, WHILE THE (ALTITUDE)
!! TOP CAN BE NEUMANN OR DIRICHLET.  BOTTOM (ALTITUDE)
!! IS ALWAYS ASSUMED TO BE DIRICHLET.


integer, dimension(4) :: flagsdirich


! This is now set in the parent procedure
!if (flagdirich==0) then      !convert current into potential normal derivative
!  Vmaxx1alt=-1*Vmaxx1*x%h1all(lx1,1:lx2all,1:lx3all)/sig0(lx1,1:lx2all,1:lx3all)
!!  print*, 'Transforming boundary from current...'
!!  print*, shape(Vmaxx1),shape(x%h1all),shape(sig0)
!else                         !Dirichlet boundary conditions, don't change
!  Vmaxx1alt=Vmaxx1
!end if


! Convert ionospheric top boundary conditions (top current or potential) into 4 sides for solver
if (flagdirich==0) then    ! Neumann top
  if (gridflag==1) then    !inverted
    flagsdirich=[0,0,1,1]
  else
    flagsdirich=[0,0,1,1]
  end if
else                       ! Dirichlet top
  if (gridflag==1) then    ! inverted
    flagsdirich=[1,0,1,1]
  else
    flagsdirich=[0,1,1,1]
  end if
end if

! FIXME: need to choose coordinates x2,x3 based on the grid size, note BC's already swapped
if (x%lx2all==1) then
  !print*, '2D, x1-x3 solve chosen'
  potential2D_fieldresolved=elliptic2D_cart(srcterm,sig0,sigP,Vminx1,Vmaxx1,Vminx3,Vmaxx3, &
                       x%dx1,x%dx1i,x%dx3all,x%dx3iall,flagsdirich,perflag,gridflag,it)
else if (x%lx3all==1) then
  !print*, '2D, x1-x2 solve chosen'
  potential2D_fieldresolved=elliptic2D_cart(srcterm,sig0,sigP,Vminx1,Vmaxx1,Vminx3,Vmaxx3, &
                       x%dx1,x%dx1i,x%dx2all,x%dx2iall,flagsdirich,perflag,gridflag,it)
end if

end procedure potential2D_fieldresolved

end submodule potential2d
