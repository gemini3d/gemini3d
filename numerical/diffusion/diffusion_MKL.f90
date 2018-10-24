module diffusion

!------ NOTE THAT SOME OLDER BACKWARD EULER ROUTINES HAVE BEEN RETAINED IN CASE
!------ THEY ARE NEEDED FOR FUTURE WORK

use phys_consts, only: kB,lsp,gammas,mindensdiv, wp
use grid, only : curvmesh,gridflag             !sizes are not imported in case we want to do subgrid diffusion
use lapack95, only: gbsv!,gtsv !banded and tridiagonal solvers
implicit none

interface TRBDF23D
  module procedure TRBDF23D_curv
end interface TRBDF23D


contains


  subroutine diffusion_prep(isp,x,lambda,betacoeff,ns,T,A,B,C,D,E,Tn,Teinf)

    !------------------------------------------------------------
    !-------COMPUTE COEFFICIENTS IN DIFFUSION EQUATIONS AND LOAD 
    !-------UP GHOST CELLS
    !------------------------------------------------------------
    !-------Note that it is done on a per species basis
    !-------Also this is never called over the full grid

    integer, intent(in) :: isp
    type(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:), intent(in) :: lambda,betacoeff

    real(wp), dimension(-1:,-1:,-1:), intent(in) :: ns
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: T
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4), &
             intent(out) :: A,B,C,D,E

    real(wp), dimension(:,:,:), intent(in) :: Tn
    real(wp), intent(in) :: Teinf

    real(wp) :: Tn0

    integer :: lx1,lx2,lx3,ix2,ix3

    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4


    !COEFFICIENTS OF PARABOLIC EQUAITON
    A(:,:,:)=0d0
    C(:,:,:)=(gammas(isp)-1d0)/kB/max(ns(1:lx1,1:lx2,1:lx3),mindensdiv)/ &
               (x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3))
    B(:,:,:)=C(:,:,:)*betacoeff/x%h1(1:lx1,1:lx2,1:lx3)    !beta must be set to zero if not electrons!
    D=lambda*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)
    E=0d0


    !SET THE BOUNDARY CONDITIONS BASED ON GRID TYPE
    if (gridflag==0) then
      do ix3=1,lx3
        do ix2=1,lx2
          Tn0=Tn(1,ix2,ix3)
          T(0,ix2,ix3)=Tn0   
          Tn0=Tn(lx1,ix2,ix3)
          T(lx1+1,ix2,ix3)=Tn0
        end do
      end do
    else if (gridflag==1) then
      do ix3=1,lx3
        do ix2=1,lx2
          Tn0=Tn(lx1,ix2,ix3)
          T(lx1+1,ix2,ix3)=Tn0   !bottom
          T(0,ix2,ix3)=Teinf     !top
        end do
      end do
    else
      do ix3=1,lx3
        do ix2=1,lx2
          Tn0=Tn(1,ix2,ix3)
          T(0,ix2,ix3)=Tn0    !bottom
          T(lx1+1,ix2,ix3)=Teinf    !top
        end do
      end do
    end if

  end subroutine diffusion_prep


  function backEuler3D(f,A,B,C,D,E,dt,dx1,dx1i)

    !------------------------------------------------------------
    !-------SOLVE A 3D SEQUENCE OF 1D DIFFUSION PROBLEMS.
    !-------GHOST CELLS ARE ACCOMODATED AS THEY PROVIDE
    !-------A CONVENIENT MEMORY SPACE FOR BOUNDARY CONDITIONS.
    !------------------------------------------------------------

    real(wp), dimension(:,:,:), intent(in) :: A,B,C,D,E   !trimmed to grid size
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: f    !expected to include ghosts
    real(wp), intent(in) :: dt
    real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference (with ghosts)
    real(wp), dimension(:), intent(in) :: dx1i   !ith centered difference

    integer :: ix1,ix2,ix3,lx1,lx2,lx3
    real(wp),dimension(size(f,1)-4) :: fx1slice

    real(wp), dimension(-1:size(f,1)-2,-1:size(f,2)-2,-1:size(f,3)-2) :: backEuler3D

    lx1=size(f,1)-4
    lx2=size(f,2)-4
    lx3=size(f,3)-4

    do ix3=1,lx3
      do ix2=1,lx2
        fx1slice=f(1:lx1,ix2,ix3)
        fx1slice=backEuler1D(fx1slice,A(:,ix2,ix3), &
                      B(:,ix2,ix3),C(:,ix2,ix3),D(:,ix2,ix3),E(:,ix2,ix3), &
                      f(0,ix2,ix3),f(lx1+1,ix2,ix3),dt,dx1,dx1i)    !inner ghost cells include boundary conditions
        backEuler3D(1:lx1,ix2,ix3)=fx1slice
      end do
    end do

  end function backEuler3D


  function TRBDF23D_curv(f,A,B,C,D,E,dt,x)

    !------------------------------------------------------------
    !-------SOLVE A 3D SEQUENCE OF 1D DIFFUSION PROBLEMS.
    !-------GHOST CELLS ARE ACCOMODATED AS THEY PROVIDE
    !-------A CONVENIENT MEMORY SPACE FOR BOUNDARY CONDITIONS.
    !------------------------------------------------------------

    real(wp), dimension(:,:,:), intent(in) :: A,B,C,D,E   !trimmed to grid size
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: f    !expected to include ghosts
    real(wp), intent(in) :: dt
    type(curvmesh), intent(in) :: x

    integer :: ix1,ix2,ix3,lx1,lx2,lx3
    real(wp),dimension(size(f,1)-4) :: fx1slice

    real(wp), dimension(-1:size(f,1)-2,-1:size(f,2)-2,-1:size(f,3)-2) :: TRBDF23D_curv

    lx1=size(f,1)-4
    lx2=size(f,2)-4
    lx3=size(f,3)-4

    do ix3=1,lx3
      do ix2=1,lx2
        fx1slice=f(1:lx1,ix2,ix3)
        fx1slice=TRBDF21D_curv(fx1slice,A(:,ix2,ix3), &
                      B(:,ix2,ix3),C(:,ix2,ix3),D(:,ix2,ix3),E(:,ix2,ix3), &
                      f(0,ix2,ix3),f(lx1+1,ix2,ix3),dt,x)             !inner ghost cells include boundary conditions
        TRBDF23D_curv(1:lx1,ix2,ix3)=fx1slice
      end do
    end do

  end function TRBDF23D_curv


  function TRBDF21D_curv(Ts,A,B,C,D,E,Tsminx1,Tsmaxx1,dt,x)

    !------------------------------------------------------------
    !-------SOLVE A 1D DIFFUSION PROBLEM.  IT IS EXPECTED THAT 
    !-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
    !-------THEY ARE PASSED INTO THIS ROUTINE.  FORM OF EQUATION
    !-------SOLVED IS:
    !------- dT/dt + A T + B dT/dx + C d/dx(D dT/dx) = E
    !------- NOTE: UPON FURTHER REVIEW I THINK THE FORM SOLVED IS
    !------- ACTUALLY:
    !------- dT/dt = A T + B dT/dx + C d/dx(D dT/dx) + E
    !------------------------------------------------------------

    real(wp), dimension(:), intent(in) :: A,B,C,D,E
    real(wp), dimension(:), intent(in) :: Ts
    real(wp), intent(in) :: Tsminx1, Tsmaxx1, dt
    type(curvmesh), intent(in) :: x
    integer, parameter :: ll=2                   !number of lower diagonals

    real(wp), dimension(3*ll+1,size(Ts)) :: M    !note extra rows for lapack workspace
    real(wp), dimension(size(Ts)) :: Dh
    integer :: ix1,lx1

    real(wp), dimension(size(Ts)) :: TR
    real(wp), dimension(size(Ts)) :: TRBDF21D_curv


    !ORGANIZE SIZES AND THERMAL CONDUCTIVITY
    lx1=size(Ts)
    Dh(1)=0.0
    Dh(2:lx1)=0.5*(D(1:lx1-1)+D(2:lx1))         !ith left cell wall thermal conductivity
!    TR(:)=Ts(:)/dt+E(:)                !boundaries to be overwritten later...  This is now done for each grid point in a separate statement


    !------------------------------------------------------------
    !-------TR HALF STEP:  DEFINE A MATRIX USING BANDED STORAGE
    !------------------------------------------------------------
    !MINX1 BOUNDARY (DIRICHLET)
    ix1=1
    M(ll+3,ix1)=1.0
    M(ll+2,ix1+1)=0.0
    M(ll+1,ix1+2)=0.0
    TR(ix1)=Tsminx1


    !FIRST INTERIOR GRID POINT
    ix1=2
    M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0 &            !ix1-1
               +B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))/2d0
    M(ll+3,ix1)=1.0/(dt/2d0)-A(ix1)/2d0 &                               !ix1
             +C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0 &
             +C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0
    M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0 &        !ix1+1, super-diag.
             -1*B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))/2d0
    M(ll+1,ix1+2)=0.0
    TR(ix1)=Ts(ix1)/(dt/2d0)+E(ix1) & 
      -M(ll+4,ix1-1)*Ts(ix1-1) &
      -(-A(ix1)/2d0+C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0+C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0)*Ts(ix1) &
      -M(ll+2,ix1+1)*Ts(ix1+1) &
      -M(ll+1,ix1+2)*Ts(ix1+2)


    !INTERIOR GRID POINTS
    forall (ix1=3:lx1-2)
      M(ll+5,ix1-2)=0.0                                               !ix1-2 grid point, sub-diag.
      M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0 &        !ix1-1
                 +B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))/2d0
      M(ll+3,ix1)=1.0/(dt/2d0)-A(ix1)/2d0 &                           !ix1
               +C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0 &
               +C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0
      M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0 &    !ix1+1, super-diag.
               -1*B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))/2d0
      M(ll+1,ix1+2)=0.0                                               !ix1+2 grid point
      TR(ix1)=Ts(ix1)/(dt/2d0)+E(ix1) &
        -M(ll+5,ix1-2)*Ts(ix1-2) & 
        -M(ll+4,ix1-1)*Ts(ix1-1) &
        -(-A(ix1)/2d0+C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0+C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0)*Ts(ix1) &
        -M(ll+2,ix1+1)*Ts(ix1+1) &
        -M(ll+1,ix1+2)*Ts(ix1+2)
    end forall


    !LAST INTERIOR GRID POINT
    ix1=lx1-1
    M(ll+5,ix1-2)=0.0
    M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0 &            !ix1-1
               +B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))/2d0
    M(ll+3,ix1)=1.0/(dt/2d0)-A(ix1)/2d0 &                               !ix1
             +C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0 &
             +C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0
    M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0 &        !ix1+1, super-diag.
             -1*B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))/2d0
    TR(ix1)=Ts(ix1)/(dt/2d0)+E(ix1) &
      -M(ll+5,ix1-2)*Ts(ix1-2) &
      -M(ll+4,ix1-1)*Ts(ix1-1) &
      -(-A(ix1)/2d0+C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1)/2d0+C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)/2d0)*Ts(ix1) &
      -M(ll+2,ix1+1)*Ts(ix1+1)


    !MAXX1 BOUNDARY
    ix1=lx1
    M(ll+5,ix1-2)=0.0
    M(ll+4,ix1-1)=0.0
    M(ll+3,ix1)=1.0
    TR(ix1)=Tsmaxx1


    !------------------------------------------------------------
    !-------TR HALF STEP MATRIX SOLUTION:  CALL LAPACK'S BANDED SOLVER
    !------------------------------------------------------------
    !BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
    call gbsv(M,TR,kl=2)



    !------------------------------------------------------------
    !-------BDF2 STEP:  DEFINE A MATRIX USING BANDED STORAGE
    !------------------------------------------------------------
    !MINX1 BOUNDARY (DIRICHLET)
    ix1=1
    M(ll+3,ix1)=1.0
    M(ll+2,ix1+1)=0.0
    M(ll+1,ix1+2)=0.0
    TRBDF21D_curv(ix1)=Tsminx1


    !FIRST INTERIOR GRID POINT
    ix1=2
    M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1) &            !ix1-1
               +B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))
    M(ll+3,ix1)=1.0/(dt/3d0)-A(ix1) &                               !ix1
             +C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1) &
             +C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)
    M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1) &        !ix1+1, super-diag.
             -1*B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))
    M(ll+1,ix1+2)=0.0
    TRBDF21D_curv(ix1)=E(ix1) & 
      -1d0/3d0*Ts(ix1)/(dt/3d0) &
      +4d0/3d0*TR(ix1)/(dt/3d0)


    !INTERIOR GRID POINTS
    forall (ix1=3:lx1-2)
      M(ll+5,ix1-2)=0.0                                               !ix1-2 grid point, sub-diag.
      M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1) &        !ix1-1
                 +B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))
      M(ll+3,ix1)=1.0/(dt/3d0)-A(ix1) &                           !ix1
               +C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1) &
               +C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)
      M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1) &    !ix1+1, super-diag.
               -1*B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))
      M(ll+1,ix1+2)=0.0                                               !ix1+2 grid point
      TRBDF21D_curv(ix1)=E(ix1) &
        -1d0/3d0*Ts(ix1)/(dt/3d0) &
        +4d0/3d0*TR(ix1)/(dt/3d0)
    end forall


    !LAST INTERIOR GRID POINT
    ix1=lx1-1
    M(ll+5,ix1-2)=0.0
    M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1) &            !ix1-1
               +B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))
    M(ll+3,ix1)=1.0/(dt/3d0)-A(ix1) &                               !ix1
             +C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1) &
             +C(ix1)*Dh(ix1)/x%dx1i(ix1)/x%dx1(ix1)
    M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/x%dx1i(ix1)/x%dx1(ix1+1) &        !ix1+1, super-diag.
             -1*B(ix1)/(x%dx1(ix1+1)+x%dx1(ix1))
    TRBDF21D_curv(ix1)=E(ix1) &
      -1d0/3d0*Ts(ix1)/(dt/3d0) &
      +4d0/3d0*TR(ix1)/(dt/3d0)

    !MAXX1 BOUNDARY
    ix1=lx1
    M(ll+5,ix1-2)=0.0
    M(ll+4,ix1-1)=0.0
    M(ll+3,ix1)=1.0
    TRBDF21D_curv(ix1)=Tsmaxx1


    !------------------------------------------------------------
    !-------BDF2 STEP MATRIX SOLUTION:  CALL LAPACK'S BANDED SOLVER
    !------------------------------------------------------------
    !BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
    call gbsv(M,TRBDF21D_curv,kl=2)

  end function TRBDF21D_curv


  function backEuler1D(Ts,A,B,C,D,E,Tsminx1,Tsmaxx1,dt,dx1,dx1i)

    !------------------------------------------------------------
    !-------SOLVE A 1D DIFFUSION PROBLEM.  IT IS EXPECTED THAT 
    !-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
    !-------THEY ARE PASSED INTO THIS ROUTINE.  FORM OF EQUATION
    !-------SOLVED IS:
    !------- dT/dt + A T + B dT/dx + C d/dx(D dT/dx) = E
    !------------------------------------------------------------

    real(wp), dimension(:), intent(in) :: A,B,C,D,E
    real(wp), dimension(:), intent(in) :: Ts
    real(wp), intent(in) :: Tsminx1, Tsmaxx1, dt
    real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
    real(wp), dimension(:), intent(in) :: dx1i   !ith centered difference
    integer, parameter :: ll=2                   !number of lower diagonals

    real(wp), dimension(3*ll+1,size(Ts)) :: M    !note extra rows for lapack workspace
    real(wp), dimension(size(Ts)) :: Dh
    integer :: ix1,lx1

    real(wp), dimension(size(Ts)) :: backEuler1D

    !------------------------------------------------------------
    !-------DEFINE A MATRIX USING BANDED STORAGE
    !------------------------------------------------------------
    lx1=size(Ts)
    Dh(1)=0.0
    Dh(2:lx1)=0.5*(D(1:lx1-1)+D(2:lx1))         !ith left cell wall thermal conductivity
    backEuler1D(:)=Ts(:)/dt+E(:)                !boundaries to be overwritten later...


    !MINX1 BOUNDARY
    ix1=1
    M(ll+3,ix1)=1.0
    M(ll+2,ix1+1)=0.0
    M(ll+1,ix1+2)=0.0
    backEuler1D(ix1)=Tsminx1


    !FIRST INTERIOR GRID POINT
    ix1=2
    M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
               +B(ix1)/(dx1(ix1+1)+dx1(ix1))
    M(ll+3,ix1)=1.0/dt-A(ix1) &                                     !ix1
             +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
             +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
    M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
             -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))
    M(ll+1,ix1+2)=0.0


    !INTERIOR GRID POINTS
    forall (ix1=3:lx1-2)
      M(ll+5,ix1-2)=0.0                                               !ix1-2 grid point, sub-diag.
      M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
                 +B(ix1)/(dx1(ix1+1)+dx1(ix1))
      M(ll+3,ix1)=1.0/dt-A(ix1) &                                     !ix1
               +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
               +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
      M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
               -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))
      M(ll+1,ix1+2)=0.0                                               !ix1+2 grid point
    end forall


    !LAST INTERIOR GRID POINT
    ix1=lx1-1
    M(ll+5,ix1-2)=0.0
    M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
               +B(ix1)/(dx1(ix1+1)+dx1(ix1))
    M(ll+3,ix1)=1.0/dt-A(ix1) &                                     !ix1
             +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
             +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
    M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
             -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))


    !MAXX1 BOUNDARY
    ix1=lx1
    M(ll+5,ix1-2)=0.0
    M(ll+4,ix1-1)=0.0
    M(ll+3,ix1)=1.0
    backEuler1D(ix1)=Tsmaxx1


    !------------------------------------------------------------
    !-------DO SOME STUFF TO CALL LAPACK'S BANDED SOLVER
    !------------------------------------------------------------
    !BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
    call gbsv(M,backEuler1D,kl=2)

  end function backEuler1D


end module diffusion


