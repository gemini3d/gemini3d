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

use phys_consts, only: wp, debug
use calculus, only: grad3D1, grad3D2, grad3D3
use meshobj, only: curvmesh
use interpolation, only: interp1
use PDEelliptic, only: elliptic3D_cart, elliptic3D_cart_periodic

implicit none (type, external)
private
public :: potential3D_fieldresolved_decimate, potential2D_polarization, potential2D_polarization_periodic, &
            potential2D_fieldresolved, potential3D_fieldresolved, potential3D_fieldresolved_truncate, &
            mumps_perm
integer, dimension(:), pointer, protected :: mumps_perm

interface ! potential2d.f90
  module function potential2D_polarization(srcterm,SigP2,SigP3,SigH,Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x,Phi0,perflag,it)
    real(wp), dimension(:,:), intent(in) :: srcterm,SigP2,SigP3,SigH,Cm,v2,v3
    !! ZZZ - THESE WILL NEED TO BE MODIFIED CONDUCTIVITIES, AND WE'LL NEED THREE OF THEM
    real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:), intent(in) :: Vminx3,Vmaxx3
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:), intent(in) :: Phi0
    logical, intent(in) :: perflag
    integer, intent(in) :: it
    real(wp), dimension(size(SigP2,1),size(SigP2,2)) :: potential2D_polarization
  end function potential2D_polarization

  module function potential2D_polarization_periodic(srcterm,SigP,SigH,Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x,Phi0,perflag,it)
    real(wp), dimension(:,:), intent(in) :: srcterm,SigP,SigH,Cm,v2,v3
    real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:), intent(in) :: Vminx3,Vmaxx3
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:), intent(in) :: Phi0
    logical, intent(in) :: perflag
    integer, intent(in) :: it
    real(wp), dimension(size(SigP,1),size(SigP,2)) :: potential2D_polarization_periodic
  end function potential2D_polarization_periodic

  module function potential2D_fieldresolved(srcterm,sig0,sigP,Vminx1,Vmaxx1,Vminx3,Vmaxx3,x,flagdirich,perflag,it)
    real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP   !arrays passed in will still have full rank 3
    real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: flagdirich
    logical, intent(in) :: perflag
    integer, intent(in) :: it
    real(wp), dimension(size(sig0,1),size(sig0,2),size(sig0,3)) :: potential2D_fieldresolved
  end function potential2D_fieldresolved
end interface

contains
  function potential3D_fieldresolved_decimate(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                    x,flagdirich,perflag,it)
    !! SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
    !! ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
    !! LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
    !! GRID POINTS.

    real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP,sigH
    real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: flagdirich
    logical, intent(in) :: perflag
    integer, intent(in) :: it

    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigP2,gradsigP3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigH2,gradsigH3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsig01
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: Ac,Bc,Cc,Dc,Ec,Fc

    integer :: lx1,lx2,lx3,ix2,ix3
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


  function potential3D_fieldresolved(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                    x,flagdirich,perflag,it)
    !! SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
    !! ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
    !! LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
    !! GRID POINTS.
    real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP,sigH
    real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: flagdirich
    logical, intent(in) :: perflag
    integer, intent(in) :: it
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigP2,gradsigP3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigH2,gradsigH3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsig01
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: Ac,Bc,Cc,Dc,Ec,Fc
    integer :: lx1,lx2,lx3
    integer, parameter :: ldec=11
    real(wp), dimension(1:size(Vminx1,1),1:size(Vminx1,2)) :: Vminx1pot,Vmaxx1pot
    real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: potential3D_fieldresolved

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

    ! Need to correct the derivatives if periodic in x3 chosen; FIXME: assume cartesian for now
    !!  FIXME: ideally this needs to be a routine from within the gradent submodule of calculus
    if (x%flagper) then
      if (debug) print*, 'Adjusting conductivity derivatives to account for periodic grid...'
      gradsigP3(:,:,1)=(sigP(:,:,2)-sigP(:,:,lx3))/(x%dx3all(2)+x%dx3all(1))
      gradsigH3(:,:,1)=(sigH(:,:,2)-sigH(:,:,lx3))/(x%dx3all(2)+x%dx3all(1))
      gradsigP3(:,:,lx3)=(sigP(:,:,1)-sigP(:,:,lx3-1))/(x%dx3all(1)+x%dx3all(lx3))
      gradsigH3(:,:,lx3)=(sigH(:,:,1)-sigH(:,:,lx3-1))/(x%dx3all(1)+x%dx3all(lx3))
    end if

    ! coefficients for 3D solve
    Ac=sigP
    Bc=sigP
    Cc=sig0
    Dc=gradsigP2+gradsigH3
    Ec=gradsigP3-gradsigH2
    Fc=gradsig01

    !ADJUST THE BOUNDARY CONDITION TO POTENTIAL DERIVATIVE INSTEAD OF CURRENT DENSITY
    Vminx1pot=-1*Vminx1/sig0(1,:,:)
    Vmaxx1pot=-1*Vmaxx1/sig0(lx1,:,:)

    ! call solver on the full grid
    if (x%flagper) then
      if (debug) print*, 'Calling 3D solve on full grid, periodic'
      potential3D_fieldresolved=elliptic3D_cart_periodic(srcterm,Ac,Bc,Cc,Dc,Ec,Fc,Vminx1pot,Vmaxx1pot, &
                      Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                      x%dx1,x%dx1i,x%dx2all,x%dx2iall,x%dx3all,x%dx3iall,flagdirich,perflag,it)
    else
      if (debug) print*, 'Calling 3D solve on full grid, non-periodic'
      potential3D_fieldresolved=elliptic3D_cart(srcterm,Ac,Bc,Cc,Dc,Ec,Fc,Vminx1pot,Vmaxx1pot, &
                      Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                      x%dx1,x%dx1i,x%dx2all,x%dx2iall,x%dx3all,x%dx3iall,flagdirich,perflag,it)
    end if
  end function potential3D_fieldresolved


  function potential3D_fieldresolved_truncate(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                    x,flagdirich,perflag,it)
    !! SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
    !! ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
    !! LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
    !! GRID POINTS.
    real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP,sigH
    real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: flagdirich
    logical, intent(in) :: perflag
    integer, intent(in) :: it
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigP2,gradsigP3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigH2,gradsigH3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsig01
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: Ac,Bc,Cc,Dc,Ec,Fc
    integer :: lx1,lx2,lx3,ix1
    real(wp), dimension(:), allocatable :: dx1trunc
    real(wp), dimension(:), allocatable :: dx1itrunc
    real(wp), dimension(:,:,:), allocatable :: Actrunc,Bctrunc,Cctrunc,Dctrunc,Ectrunc,Fctrunc,srctermtrunc
    real(wp), dimension(:,:), allocatable :: Vminx2trunc,Vmaxx2trunc
    real(wp), dimension(:,:), allocatable :: Vminx3trunc, Vmaxx3trunc
    real(wp), dimension(:,:,:), allocatable :: Phitrunc
    real(wp), dimension(1:size(Vminx1,1),1:size(Vminx1,2)) :: Vminx1pot,Vmaxx1pot
    real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: potential3D_fieldresolved_truncate
    integer :: lx1trunc
    real(wp) :: alttrunc

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

    ! Need to correct the derivatives if periodic in x3 chosen; FIXME: assume cartesian for now
    !!  FIXME: ideally this needs to be a routine from within the gradent submodule of calculus
    if (x%flagper) then
      if (debug) print*, 'Adjusting conductivity derivatives to account for periodic grid...'
      gradsigP3(:,:,1)=(sigP(:,:,2)-sigP(:,:,lx3))/(x%dx3all(2)+x%dx3all(1))
      gradsigH3(:,:,1)=(sigH(:,:,2)-sigH(:,:,lx3))/(x%dx3all(2)+x%dx3all(1))
      gradsigP3(:,:,lx3)=(sigP(:,:,1)-sigP(:,:,lx3-1))/(x%dx3all(1)+x%dx3all(lx3))
      gradsigH3(:,:,lx3)=(sigH(:,:,1)-sigH(:,:,lx3-1))/(x%dx3all(1)+x%dx3all(lx3))
    end if

    ! coefficients for 3D solve
    Ac=sigP
    Bc=sigP
    Cc=sig0
    Dc=gradsigP2+gradsigH3
    Ec=gradsigP3-gradsigH2
    Fc=gradsig01

    ! Now truncate the problem at some appropriate altitude
    ix1=1
    alttrunc=300e3
    do while(x%alt(ix1,1,1)<alttrunc)
      ix1=ix1+1
    end do
    if (ix1>=1) then
      lx1trunc=ix1-1
      print*, 'Truncating field aligned solve at index:  ',lx1trunc
    else
      print*, 'Unable to truncate field-resolved solve at sensible altitude...',ix1,alttrunc
      error stop
    end if

    !ADJUST THE BOUNDARY CONDITION TO POTENTIAL DERIVATIVE INSTEAD OF CURRENT DENSITY
    Vminx1pot=-1*Vminx1/sig0(1,:,:)
    Vmaxx1pot=-1*Vmaxx1/sig0(lx1trunc,:,:)

    ! allocate and assign truncated variables
    allocate(Actrunc(1:lx1trunc,1:lx2,1:lx3))
    allocate(Bctrunc,Cctrunc,Dctrunc,Ectrunc,Fctrunc,srctermtrunc, mold=Actrunc)
    allocate(Vminx2trunc(1:lx1trunc,1:lx3))
    allocate(Vmaxx2trunc, mold=Vminx2trunc)
    allocate(Vminx3trunc(1:lx1trunc,1:lx2))
    allocate(Vmaxx3trunc, mold=Vminx3trunc)
    allocate(dx1trunc(0:lx1trunc+2))
    allocate(dx1itrunc(1:lx1trunc+1))
    Actrunc=Ac(1:lx1trunc,1:lx2,1:lx3)
    Bctrunc=Bc(1:lx1trunc,1:lx2,1:lx3)
    Cctrunc=Cc(1:lx1trunc,1:lx2,1:lx3)
    Dctrunc=Dc(1:lx1trunc,1:lx2,1:lx3)
    Ectrunc=Ec(1:lx1trunc,1:lx2,1:lx3)
    Fctrunc=Fc(1:lx1trunc,1:lx2,1:lx3)
    srctermtrunc=srcterm(1:lx1trunc,1:lx2,1:lx3)
    Vminx2trunc=Vminx2(1:lx1trunc,1:lx3)
    Vmaxx2trunc=Vmaxx2(1:lx1trunc,1:lx3)
    Vminx3trunc=Vminx3(1:lx1trunc,1:lx2)
    Vmaxx3trunc=Vmaxx3(1:lx1trunc,1:lx2)
    dx1trunc=x%dx1(0:lx1trunc+2)
    dx1itrunc=x%dx1i(1:lx1trunc+1)

    ! call the solver for the truncated grid
    allocate(Phitrunc, mold=Actrunc)
    if (x%flagper) then
      if (debug) print*, 'Calling 3D solve on full grid, periodic'
      Phitrunc=elliptic3D_cart_periodic(srctermtrunc,Actrunc,Bctrunc,Cctrunc,Dctrunc,Ectrunc,Fctrunc,Vminx1pot,Vmaxx1pot, &
                      Vminx2trunc,Vmaxx2trunc,Vminx3trunc,Vmaxx3trunc, &
                      dx1trunc,dx1itrunc,x%dx2all,x%dx2iall,x%dx3all,x%dx3iall,flagdirich,perflag,it)
    else
      if (debug) print*, 'Calling 3D solve on full grid, non-periodic'
      Phitrunc=elliptic3D_cart(srctermtrunc,Actrunc,Bctrunc,Cctrunc,Dctrunc,Ectrunc,Fctrunc,Vminx1pot,Vmaxx1pot, &
                      Vminx2trunc,Vmaxx2trunc,Vminx3trunc,Vmaxx3trunc, &
                      dx1trunc,dx1itrunc,x%dx2all,x%dx2iall,x%dx3all,x%dx3iall,flagdirich,perflag,it)
    end if
    potential3D_fieldresolved_truncate(1:lx1trunc,1:lx2,1:lx3)=Phitrunc(1:lx1trunc,1:lx2,1:lx3)
    do ix1=lx1trunc+1,lx1
      potential3D_fieldresolved_truncate(ix1,1:lx2,1:lx3)=Phitrunc(lx1trunc,1:lx2,1:lx3)
    end do

    ! deallocate temporary arrays
    deallocate(Actrunc,Bctrunc,Cctrunc,Dctrunc,Ectrunc,Fctrunc,srctermtrunc,Vminx2trunc,Vmaxx2trunc, &
               Vminx3trunc,Vmaxx3trunc,dx1trunc,dx1itrunc,Phitrunc)
  end function potential3D_fieldresolved_truncate
end module potential_mumps
