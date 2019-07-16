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

use phys_consts, only: wp
use calculus, only: grad3D1, grad3D2, grad3D3, grad2D1_curv_alt, grad2D3, grad2D3_curv_periodic
use grid, only: curvmesh, gridflag
use mpimod, only: myid
use interpolation, only:interp1

implicit none

private

#if REALBITS==32
include 'smumps_struc.h'
#elif REALBITS==64
include 'dmumps_struc.h'
#else
error stop "realbits must be 32 or 64"
#endif

integer, dimension(:), pointer, protected, save :: mumps_perm   !cached permutation, unclear whether save is necessary...

public :: elliptic3d_curv, elliptic3D_decimate, elliptic2d_pol_conv_curv, elliptic2d_pol_conv_curv_periodic2, &
          elliptic2d_nonint_curv, elliptic_workers

contains


function elliptic3D_decimate(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
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
integer :: ldec
real(wp), dimension(:), allocatable :: x1dec
real(wp), dimension(:), allocatable :: dx1dec
real(wp), dimension(:), allocatable :: x1idec
real(wp), dimension(:), allocatable :: dx1idec
real(wp), dimension(:,:,:), allocatable :: Acdec,Bcdec,Ccdec,Dcdec,Ecdec,Fcdec,srctermdec
real(wp), dimension(:,:), allocatable :: Vminx2dec,Vmaxx2dec
real(wp), dimension(:,:), allocatable :: Vminx3dec, Vmaxx3dec
real(wp), dimension(:,:,:), allocatable :: Phidec

real(wp), dimension(1:size(Vminx1,1),1:size(Vminx1,2)) :: Vminx1pot,Vmaxx1pot

real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D_decimate



!SYSTEM SIZES
lx1=size(srcterm,1)    !These will be full grid sizes if called from root (only acceptable thing)
lx2=size(srcterm,2)
lx3=size(srcterm,3)


!COMPUTE AUXILIARY COEFFICIENTS TO PASS TO CART SOLVER
print *, 'Prepping coefficients for elliptic equation...'
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
print*, 'Decimating parallel grid...'
ldec=11
allocate(x1dec(-1:ldec+2),dx1dec(0:ldec+2),x1idec(ldec+1),dx1idec(ldec))
x1dec(-1:lx1+2)=[x%x1(-1),x%x1(0),x%x1(1),81.8e3_wp,84.2e3_wp,87.5e3_wp,93.3e3_wp,106.0e3_wp,124.0e3_wp, &
            144.6e3_wp,206.7e3_wp,882.2e3_wp,x%x1(lx1),x%x1(lx1+1),x%x1(lx1+2)]
dx1dec(0:ldec+2)=x1dec(0:ldec+2)-x1dec(-1:ldec+1)
x1idec(1:ldec+1)=0.5_wp*(x1dec(0:ldec)+x1dec(1:ldec+1))
dx1idec(1:ldec)=x1idec(2:ldec+1)-x1idec(1:ldec)


!INTERPOLATE COEFFICIENTS AND SOURCE TERM ONTO DECIMATED GRID
print*, 'Interpolating coefficients...'
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
print*, 'Calling solve on decimated grid...'
allocate(Phidec(1:ldec,1:lx2,1:lx3))
Phidec=elliptic3D_cart(srctermdec,Acdec,Bcdec,Ccdec,Dcdec,Ecdec,Fcdec,Vminx1pot,Vmaxx1pot, &
                Vminx2dec,Vmaxx2dec,Vminx3dec,Vmaxx3dec, &
                dx1dec,dx1idec,x%dx2all,x%dx2iall,x%dx3all,x%dx3iall,flagdirich,perflag,it)


!INTERPOLATE BACK UP TO MAIN GRID
print*, 'Upsampling potential...'
do ix2=1,lx2
  do ix3=1,lx3
    elliptic3D_decimate(:,ix2,ix3)=interp1(x1dec,Phidec(:,ix2,ix3),x%x1(1:lx1))
  end do
end do


!AGAIN NEED TO FIX THE EDGES...
elliptic3D_decimate(1,:,:)=Phidec(1,:,:)
elliptic3D_decimate(lx1,:,:)=Phidec(ldec,:,:)


!CLEAN UP THE ALLOCATED ARRAYS
deallocate(srctermdec,Acdec,Bcdec,Ccdec,Dcdec,Ecdec,Fcdec,dx1dec,x1dec,x1idec,dx1idec,Vminx2dec,Vmaxx2dec,Vminx3dec,Vmaxx3dec)
deallocate(Phidec)

end function elliptic3D_decimate


function elliptic3D_cart(srcterm,Ac,Bc,Cc,Dc,Ec,Fc,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                  dx1,dx1i,dx2all,dx2iall,dx3all,dx3iall,flagdirich,perflag,it)

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
!-------ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
!-------LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
!-------GRID POINTS.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: srcterm,Ac,Bc,Cc,Dc,Ec,Fc
real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
real(wp), dimension(0:), intent(in) :: dx1         !backweard diffs start at index zero due to ghost cells
real(wp), dimension(:), intent(in) :: dx1i         !centered diffs do not include any ghost cells
real(wp), dimension(0:), intent(in) :: dx2all
real(wp), dimension(:), intent(in) :: dx2iall
real(wp), dimension(0:), intent(in) :: dx3all
real(wp), dimension(:), intent(in) :: dx3iall
integer, intent(in) :: flagdirich
logical, intent(in) :: perflag
integer, intent(in) :: it

integer :: ix1,ix2,ix3,lx1,lx2,lx3
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
real(wp) :: tstart,tfin

#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#elif REALBITS==64
type (DMUMPS_STRUC) mumps_par
#else
error stop "realbits must be 32 or 64"
#endif

real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D_cart


!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
if (myid==0) then
  lx1=size(Ac,1)
  lx2=size(Ac,2)
  lx3=size(Ac,3)
  lPhi=lx1*lx2*lx3

  lent=7*(lx1-2)*(lx2-2)*(lx3-2)                                                 !interior entries
  lent=lent+2*(lx1-2)*(lx2-2)+2*(lx2-2)*(lx3-2)+2*(lx1-2)*(lx3-2)                !6 faces of cube
  lent=lent+4*(lx1-2)+4*(lx2-2)+4*(lx3-2)                                        !12 edges
  lent=lent+8                                                                    !8 corners, now we have total nonzero entries
  lent=lent+lx2*lx3                                                              !entries to deal with Neumann conditions on bottom
  if (flagdirich==0) then                                                        !more entries if Neumann on top
    lent=lent+lx2*lx3
  end if

  allocate(ir(lent),ic(lent),M(lent),b(lPhi))

  print *, 'MUMPS will attempt a solve of size:  ',lx1,lx2,lx3
  print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent



  !------------------------------------------------------------
  !-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
  !-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
  !-------GUIDE).
  !------------------------------------------------------------
  !LOAD UP MATRIX ELEMENTS
  M(:)=0d0
  b=pack(srcterm,.true.)           !boundaries overwritten later
  ient=1
  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
        iPhi=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix3) as a column vector.  Also row of big matrix

        if (ix1==1) then          !BOTTOM GRID POINTS + CORNER, USE NEUMANN HERE, PRESUMABLY ZERO
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=-1d0
          ient=ient+1
          ir(ient)=iPhi
          ic(ient)=iPhi+1
          M(ient)=1d0
  !            b(iPhi)=Vminx1(ix3)
          b(iPhi)=0d0    !force bottom current to zero
          ient=ient+1
        elseif (ix1==lx1) then    !TOP GRID POINTS + CORNER
          if (flagdirich/=0) then
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1d0
            b(iPhi)=Vmaxx1(ix2,ix3)
            ient=ient+1
          else
            ir(ient)=iPhi
            ic(ient)=iPhi-1
            M(ient)=-1d0/dx1(lx1)
            ient=ient+1
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1d0/dx1(lx1)
            b(iPhi)=Vmaxx1(ix2,ix3)
            ient=ient+1
          end if
        elseif (ix2==1) then      !LEFT BOUNDARY
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1.0
          b(iPhi)=Vminx2(ix1,ix3)
          ient=ient+1
        elseif (ix2==lx2) then    !RIGHT BOUNDARY
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1.0
          b(iPhi)=Vmaxx2(ix1,ix3)
          ient=ient+1
        elseif (ix3==1) then
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1.0
          b(iPhi)=Vminx3(ix1,ix2)
          ient=ient+1
        elseif (ix3==lx3) then
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1.0
          b(iPhi)=Vmaxx3(ix1,ix2)
          ient=ient+1
        else                      !INTERIOR
          !ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
          ir(ient)=iPhi
          ic(ient)=iPhi-lx1*lx2
          M(ient)=Bc(ix1,ix2,ix3)/dx3all(ix3)/dx3iall(ix3)-Ec(ix1,ix2,ix3)/(dx3all(ix3+1)+dx3all(ix3))
          ient=ient+1

          !ix1,ix2-1,ix3
          ir(ient)=iPhi
          ic(ient)=iPhi-lx1
          M(ient)=Ac(ix1,ix2,ix3)/dx2all(ix2)/dx2iall(ix2)-Dc(ix1,ix2,ix3)/(dx2all(ix2+1)+dx2all(ix2))
          ient=ient+1

          !ix1-1,ix2,ix3
          ir(ient)=iPhi
          ic(ient)=iPhi-1
          M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1)/dx1i(ix1)-Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))
          ient=ient+1

          !ix1,ix2,ix3
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=-1d0*Ac(ix1,ix2,ix3)*(1d0/dx2all(ix2+1)/dx2iall(ix2)+1d0/dx2all(ix2)/dx2iall(ix2))- &
                       Bc(ix1,ix2,ix3)*(1d0/dx3all(ix3+1)/dx3iall(ix3)+1d0/dx3all(ix3)/dx3iall(ix3))- &
                       Cc(ix1,ix2,ix3)*(1d0/dx1(ix1+1)/dx1i(ix1)+1d0/dx1(ix1)/dx1i(ix1))
          ient=ient+1

          !ix1+1,ix2,ix3
          ir(ient)=iPhi
          ic(ient)=iPhi+1
          M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1+1)/dx1i(ix1)+Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))
          ient=ient+1

          !ix1,ix2+1,ix3
          ir(ient)=iPhi
          ic(ient)=iPhi+lx1
          M(ient)=Ac(ix1,ix2,ix3)/dx2all(ix2+1)/dx2iall(ix2)+Dc(ix1,ix2,ix3)/(dx2all(ix2+1)+dx2all(ix2))
          ient=ient+1

          !ix1,ix2,ix3+1
          ir(ient)=iPhi
          ic(ient)=iPhi+lx1*lx2
          M(ient)=Bc(ix1,ix2,ix3)/dx3all(ix3+1)/dx3iall(ix3)+Ec(ix1,ix2,ix3)/(dx3all(ix3+1)+dx3all(ix3))
          ient=ient+1
        end if
      end do
    end do
  end do
end if
print *, 'Number of entries used:  ',ient-1


! INIT MUMPS
! MUMPS 4.10 disregarded ICNTL(4)
mumps_par%ICNTL(1) = stderr   ! error msg stream
mumps_par%ICNTL(2) = stdout   ! stats and warning stream
mumps_par%ICNTL(3) = 0  ! global information verbosity  0 = off
mumps_par%ICNTL(4) = 1  ! 1: error only   2: errors, warnings, stats.

mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif



!LOAD OUR PROBLEM
if ( myid==0 ) then
  mumps_par%N=lPhi
  mumps_par%NZ=lent
  allocate( mumps_par%IRN ( mumps_par%NZ ) )
  allocate( mumps_par%JCN ( mumps_par%NZ ) )
  allocate( mumps_par%A( mumps_par%NZ ) )
  allocate( mumps_par%RHS ( mumps_par%N  ) )
  mumps_par%IRN=ir
  mumps_par%JCN=ic
  mumps_par%A=M
  mumps_par%RHS=b
  deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

  if (perflag .and. it/=1) then       !used cached permutation
    allocate(mumps_par%PERM_IN(mumps_par%N))
    mumps_par%PERM_IN=mumps_perm
    mumps_par%ICNTL(7)=1
  end if


  !3D solves very often need better memory relaxation
  mumps_par%ICNTL(14)=200
end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
if ( myid==0 ) then
  print *, 'Now organizing results...'

  if (perflag .and. it==1) then
    allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
    mumps_perm=mumps_par%SYM_PERM
  end if

  elliptic3D_cart=reshape(mumps_par%RHS,[lx1,lx2,lx3])

  print *, 'Now attempting deallocations...'

  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
end if

mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif

end function elliptic3D_cart


function elliptic3D_curv(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
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

integer :: ix1,ix2,ix3,lx1,lx2,lx3
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
real(wp) :: tstart,tfin

#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#elif REALBITS==64
type (DMUMPS_STRUC) mumps_par
#else
error stop "realbits must be 32 or 64"
#endif

real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D_curv

!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
if (myid==0) then
  lx1=size(sig0,1)
  lx2=size(sig0,2)
  lx3=size(sig0,3)
  lPhi=lx1*lx2*lx3

  lent=7*(lx1-2)*(lx2-2)*(lx3-2)                                                 !interior entries
  lent=lent+2*(lx1-2)*(lx2-2)+2*(lx2-2)*(lx3-2)+2*(lx1-2)*(lx3-2)                !6 faces of cube
  lent=lent+4*(lx1-2)+4*(lx2-2)+4*(lx3-2)                                        !12 edges
  lent=lent+8                                                                    !8 corners, now we have total nonzero entries
  lent=lent+lx2*lx3                                                              !entries to deal with Neumann conditions on bottom
  if (flagdirich==0) then                                                        !more entries if Neumann on top
    lent=lent+lx2*lx3
  end if

  allocate(ir(lent),ic(lent),M(lent),b(lPhi))

  print *, 'MUMPS will attempt a solve of size:  ',lx1,lx2,lx3
  print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


  !COMPUTE AUXILIARY COEFFICIENTS
  print *, 'Prepping coefficients for elliptic equation...'
  gradsig01=grad3D1(sig0,x,1,lx1,1,lx2,1,lx3)
  gradsigP2=grad3D2(sigP,x,1,lx1,1,lx2,1,lx3)
  gradsigP3=grad3D3(sigP,x,1,lx1,1,lx2,1,lx3)
  gradsigH2=grad3D2(sigH,x,1,lx1,1,lx2,1,lx3)
  gradsigH3=grad3D3(sigH,x,1,lx1,1,lx2,1,lx3)

  Ac=sigP; Bc=sigP; Cc=sig0;
  Dc=gradsigP2+gradsigH3
  Ec=gradsigP3-gradsigH2
  Fc=gradsig01


  !------------------------------------------------------------
  !-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
  !-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
  !-------GUIDE).
  !------------------------------------------------------------
  !LOAD UP MATRIX ELEMENTS
  M(:)=0d0
  b=pack(srcterm,.true.)           !boundaries overwritten later
  ient=1
  do ix3=1,lx3
    do ix2=1,lx2
    do ix1=1,lx1
      iPhi=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix3) as a column vector.  Also row of big matrix

      if (ix1==1) then          !BOTTOM GRID POINTS + CORNER, USE NEUMANN HERE, PRESUMABLY ZERO
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-1d0
        ient=ient+1
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient)=1d0
!            b(iPhi)=Vminx1(ix3)
        b(iPhi)=0d0    !force bottom current to zero
        ient=ient+1
      elseif (ix1==lx1) then    !TOP GRID POINTS + CORNER
        if (flagdirich/=0) then
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1d0
          b(iPhi)=Vmaxx1(ix2,ix3)
          ient=ient+1
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-1
          M(ient)=-1d0/x%dx1(lx1)
          ient=ient+1
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1d0/x%dx1(lx1)
          b(iPhi)=Vmaxx1(ix2,ix3)
          ient=ient+1
        end if
      elseif (ix2==1) then      !LEFT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vminx2(ix1,ix3)
        ient=ient+1
      elseif (ix2==lx2) then    !RIGHT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vmaxx2(ix1,ix3)
        ient=ient+1
      elseif (ix3==1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vminx3(ix1,ix2)
        ient=ient+1
      elseif (ix3==lx3) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vmaxx3(ix1,ix2)
        ient=ient+1
      else                      !INTERIOR
        !ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1*lx2
        M(ient)=Bc(ix1,ix2,ix3)/x%dx3all(ix3)/x%dx3iall(ix3)-Ec(ix1,ix2,ix3)/(x%dx3all(ix3+1)+x%dx3all(ix3))
        ient=ient+1

        !ix1,ix2-1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1
        M(ient)=Ac(ix1,ix2,ix3)/x%dx2(ix2)/x%dx2i(ix2)-Dc(ix1,ix2,ix3)/(x%dx2(ix2+1)+x%dx2(ix2))
        ient=ient+1

        !ix1-1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-1
        M(ient)=Cc(ix1,ix2,ix3)/x%dx1(ix1)/x%dx1i(ix1)-Fc(ix1,ix2,ix3)/(x%dx1(ix1+1)+x%dx1(ix1))
        ient=ient+1

        !ix1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-1d0*Ac(ix1,ix2,ix3)*(1d0/x%dx2(ix2+1)/x%dx2i(ix2)+1d0/x%dx2(ix2)/x%dx2i(ix2))- &
                     Bc(ix1,ix2,ix3)*(1d0/x%dx3all(ix3+1)/x%dx3iall(ix3)+1d0/x%dx3all(ix3)/x%dx3iall(ix3))- &
                     Cc(ix1,ix2,ix3)*(1d0/x%dx1(ix1+1)/x%dx1i(ix1)+1d0/x%dx1(ix1)/x%dx1i(ix1))
        ient=ient+1

        !ix1+1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient)=Cc(ix1,ix2,ix3)/x%dx1(ix1+1)/x%dx1i(ix1)+Fc(ix1,ix2,ix3)/(x%dx1(ix1+1)+x%dx1(ix1))
        ient=ient+1

        !ix1,ix2+1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1
        M(ient)=Ac(ix1,ix2,ix3)/x%dx2(ix2+1)/x%dx2i(ix2)+Dc(ix1,ix2,ix3)/(x%dx2(ix2+1)+x%dx2(ix2))
        ient=ient+1

        !ix1,ix2,ix3+1
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1*lx2
        M(ient)=Bc(ix1,ix2,ix3)/x%dx3all(ix3+1)/x%dx3iall(ix3)+Ec(ix1,ix2,ix3)/(x%dx3all(ix3+1)+x%dx3all(ix3))
        ient=ient+1
      end if
    end do
  end do
  end do
end if
print *, 'Number of entries used:  ',ient-1


! INIT MUMPS
! MUMPS 4.10 disregarded ICNTL(4)
mumps_par%ICNTL(1) = stderr   ! error msg stream
mumps_par%ICNTL(2) = stdout   ! stats and warning stream
mumps_par%ICNTL(3) = 0  ! global information verbosity  0 = off
mumps_par%ICNTL(4) = 1  ! 1: error only   2: errors, warnings, stats.

mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif



!LOAD OUR PROBLEM
if ( myid==0 ) then
  mumps_par%N=lPhi
  mumps_par%NZ=lent
  allocate( mumps_par%IRN ( mumps_par%NZ ) )
  allocate( mumps_par%JCN ( mumps_par%NZ ) )
  allocate( mumps_par%A( mumps_par%NZ ) )
  allocate( mumps_par%RHS ( mumps_par%N  ) )
  mumps_par%IRN=ir
  mumps_par%JCN=ic
  mumps_par%A=M
  mumps_par%RHS=b
  deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

  if (perflag .and. it/=1) then       !used cached permutation
    allocate(mumps_par%PERM_IN(mumps_par%N))
    mumps_par%PERM_IN=mumps_perm
    mumps_par%ICNTL(7)=1
  end if
end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
if ( myid==0 ) then
  print *, 'Now organizing results...'

  if (perflag .and. it==1) then
    allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
    mumps_perm=mumps_par%SYM_PERM
  end if

  elliptic3D_curv=reshape(mumps_par%RHS,[lx1,lx2,lx3])

  print *, 'Now attempting deallocations...'

  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
end if

mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif

end function elliptic3D_curv


function elliptic2D_pol_conv_curv(srcterm,SigP2,SigP3,SigH,Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x,Phi0,perflag,it)

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

real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: SigPh2    !I'm too lazy to recode these as SigP2h2, etc.
real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: SigPh3
!    real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: SigHh2
!    real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: SigHh3
real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: Cmh2
real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: Cmh3
real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: gradSigH2,gradSigH3

real(wp) :: coeff    !coefficient for calculating polarization terms
integer :: ix2,ix3,lx2,lx3    !this overwrites the
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
real(wp) :: tstart,tfin

integer :: utrace

#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#elif REALBITS==64
type (DMUMPS_STRUC) mumps_par
#else
error stop "realbits must be 32 or 64"
#endif

real(wp), dimension(size(SigP2,1),size(SigP2,2)) :: elliptic2D_pol_conv_curv


!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
if (myid==0) then
  lx2=size(SigP2,1)    !note that these are full-grid sizes since grid module globals are not in scope
  lx3=size(SigP2,2)
  lPhi=lx2*lx3
!      lent=5*(lx2-2)*(lx3-2)+2*lx2+2*(lx3-2)    !static model; left here as a reference
  lent=17*(lx2-2)*(lx3-2)+2*lx2+2*(lx3-2)-3*2*(lx2-2)-3*2*(lx3-2)    !interior+boundary-x3_adj-x2_adj.  Note that are 3 sets of entries for each adjacent point
  allocate(ir(lent),ic(lent),M(lent),b(lPhi))

  print *, 'MUMPS will attempt a solve of size:  ',lx2,lx3
  print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


  !PREP INPUT DATA FOR SOLUTION OF SYSTEM
  SigPh2(1,:)=0d0
  SigPh2(2:lx2,:)=0.5d0*(SigP2(1:lx2-1,:)+SigP2(2:lx2,:))   !note the different conductiances here ot be associated with derivatives in different directions
  SigPh3(:,1)=0d0
  SigPh3(:,2:lx3)=0.5d0*(SigP3(:,1:lx3-1)+SigP3(:,2:lx3))
  Cmh2(1,:)=0d0
  Cmh2(2:lx2,:)=0.5d0*(Cm(1:lx2-1,:)+Cm(2:lx2,:))
  Cmh3(:,1)=0d0
  Cmh3(:,2:lx3)=0.5d0*(Cm(:,1:lx3-1)+Cm(:,2:lx3))



  !gradSigH2=grad2D1(SigH,x,1,lx2)   !x2 is now 1st index and x3 is second...  This one appears to be the problem.  This issue here is that grad2D1 automatically uses x%dx1 as the differential element...
  gradSigH2=grad2D1_curv_alt(SigH,x,1,lx2)
  !! note the alt since we need to use dx2 as differential...  Tricky bug/feature
  gradSigH3=grad2D3(SigH,x,1,lx3)
  !! awkward way of handling this special case derivative which uses x3 as the differential to operate on a 2D array.

!
!  open(newunit=utrace, form='unformatted', access='stream', file='matrixentryinfo.raw8', status='replace', action='write')
!  write(utrace) gradSigH2
!  write(utrace) gradSigH3
!  write(utrace) SigPh2
!  write(utrace) SigPh3
!  write(utrace) Cmh2
!  write(utrace) Cmh3
!  write(utrace) SigP2
!  write(utrace) SigP3
!  write(utrace) SigH
!  write(utrace) Cm
!  write(utrace) v2
!  write(utrace) v3
!  close(utrace)
!

  !------------------------------------------------------------
  !-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
  !-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
  !-------GUIDE).
  !------------------------------------------------------------

  !LOAD UP MATRIX ELEMENTS
  M(:)=0d0
  b=pack(srcterm,.true.)           !boundaries overwritten later, polarization terms also added later.
  ient=1
  do ix3=1,lx3
    do ix2=1,lx2
      iPhi=lx2*(ix3-1)+ix2     !linear index referencing Phi(ix2,ix3) as a column vector.  Also row of big matrix

      if (ix2==1) then          !BOTTOM GRID POINTS + CORNER
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1d0
        b(iPhi)=Vminx2(ix3)
        ient=ient+1
      elseif (ix2==lx2) then    !TOP GRID POINTS + CORNER
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1d0
        b(iPhi)=Vmaxx2(ix3)
        ient=ient+1
      elseif (ix3==1) then      !LEFT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1d0
        b(iPhi)=Vminx3(ix2)
        ient=ient+1
      elseif (ix3==lx3) then    !RIGHT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1d0
        b(iPhi)=Vmaxx3(ix2)
        ient=ient+1
      else                      !INTERIOR LOCATION
        !ix2-1,ix3-2 grid point
        coeff=-1d0*Cm(ix2,ix3-1)*v2(ix2,ix3-1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)+x%dx3all(ix3))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        if (ix3==2) then    !out of bounds, use nearest BC, and add to known vector
          b(iPhi)=b(iPhi)-coeff*Vminx3(ix2-1)
        else    !in bounds, add to matrix
          ir(ient)=iPhi
          ic(ient)=iPhi-2*lx2-1

          M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

          ient=ient+1
        end if


        !ix2,ix3-2 grid point
        coeff=-1d0*Cm(ix2,ix3-1)*v3(ix2,ix3-1)/( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)*x%dx3iall(ix3-1)) )
        if (ix3==2) then    !bit of intentional code duplication here and in the following sections to keep things organized in a way I can debug...
          b(iPhi)=b(iPhi)-coeff*Vminx3(ix2)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-2*lx2

          M(ient)=coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

          ient=ient+1
        end if


        !ix2+1,ix3-2 grid point
        coeff=Cm(ix2,ix3-1)*v2(ix2,ix3-1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)+x%dx3all(ix3))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        if (ix3==2) then
          b(iPhi)=b(iPhi)-coeff*Vminx3(ix2+1)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-2*lx2+1

          M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

          ient=ient+1
        end if


        !ix2-2,ix3-1
        coeff=-1d0*Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)+x%dx2all(ix2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        if (ix2==2) then
          b(iPhi)=b(iPhi)-coeff*Vminx2(ix3-1)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-lx2-2

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !!!ix2,ix3-1 grid point in ix2,ix3 equation
        ir(ient)=iPhi
        ic(ient)=iPhi-lx2

        M(ient)=SigPh3(ix2,ix3)/(x%dx3iall(ix3)*x%dx3all(ix3))+gradSigH2(ix2,ix3)/(x%dx3all(ix3)+x%dx3all(ix3+1))    !static terms

        coeff=Cmh3(ix2,ix3)/(dt*x%dx3iall(ix3)*x%dx3all(ix3))
        M(ient)=M(ient)+coeff   !polarization time derivative terms
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2,ix3-1)   !add in polarziation terms that include previous time step potential at this grid point

        coeff=Cm(ix2,ix3-1)*v3(ix2,ix3-1)/( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3)*x%dx3iall(ix3-1)) )+ &
              Cm(ix2,ix3-1)*v3(ix2,ix3-1)/( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)*x%dx3iall(ix3-1)) )
        M(ient)=M(ient)+coeff   !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

        coeff=Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)+x%dx2all(ix2+2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )+ &
              Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)+x%dx2all(ix2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        M(ient)=M(ient)+coeff   !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

        ient=ient+1


        !ix2+2,ix3-1 grid point
        coeff=-1d0*Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
              ( (x%dx2all(ix2+1)+x%dx2all(ix2+2))*(x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        if (ix2==lx2-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx2(ix3-1)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-lx2+2

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !ix2-2,ix3 grid point
        coeff=-1d0*Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)*x%dx2iall(ix2-1)) )
        if (ix2==2) then
          b(iPhi)=b(iPhi)-coeff*Vminx2(ix3)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-2

          M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

          ient=ient+1
        end if


        !!!ix2-1,ix3 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi-1

        M(ient)=SigPh2(ix2,ix3)/(x%dx2iall(ix2)*x%dx2all(ix2))-gradSigH3(ix2,ix3)/(x%dx2all(ix2)+x%dx2all(ix2+1))    !static

        coeff=Cmh2(ix2,ix3)/(dt*x%dx2iall(ix2)*x%dx2all(ix2))
        M(ient)=M(ient)+coeff    !pol. time deriv.
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2-1,ix3)    !BC's and pol. time deriv.

        coeff=Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2)*x%dx2iall(ix2-1)) )+ &
              Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)*x%dx2iall(ix2-1)) )
        M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

        coeff=Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )+ &
              Cm(ix2,ix3-1)*v2(ix2,ix3-1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)+x%dx3all(ix3))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        M(ient)=M(ient)+coeff     !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

        ient=ient+1


        !!!ix2,ix3 grid point (main diagonal)
        ir(ient)=iPhi
        ic(ient)=iPhi

        M(ient)=-1d0*SigPh2(ix2+1,ix3)/(x%dx2iall(ix2)*x%dx2all(ix2+1)) &
                -1d0*SigPh2(ix2,ix3)/(x%dx2iall(ix2)*x%dx2all(ix2)) &
                -1d0*SigPh3(ix2,ix3+1)/(x%dx3iall(ix3)*x%dx3all(ix3+1)) &
                -1d0*SigPh3(ix2,ix3)/(x%dx3iall(ix3)*x%dx3all(ix3))    !static

        coeff=-1d0*Cmh2(ix2+1,ix3)/(dt*x%dx2iall(ix2)*x%dx2all(ix2+1)) &
              -1d0*Cmh2(ix2,ix3)/(dt*x%dx2iall(ix2)*x%dx2all(ix2)) &
              -1d0*Cmh3(ix2,ix3+1)/(dt*x%dx3iall(ix3)*x%dx3all(ix3+1)) &
              -1d0*Cmh3(ix2,ix3)/(dt*x%dx3iall(ix3)*x%dx3all(ix3))
        M(ient)=M(ient)+coeff    !pol. time deriv.
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2,ix3)    !BC's and pol. time deriv.

        coeff=Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)*x%dx2iall(ix2+1)) )+ &
              (-1d0)*Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2)*x%dx2iall(ix2-1)) )
        M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

        coeff=Cm(ix2,ix3+1)*v3(ix2,ix3+1)/( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)*x%dx3all(ix3+1)) )+ &
              (-1d0)*Cm(ix2,ix3-1)*v3(ix2,ix3-1)/( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3)*x%dx3iall(ix3-1)) )
        M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

        ient=ient+1


        !!!ix2+1,ix3 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi+1

        M(ient)=SigPh2(ix2+1,ix3)/(x%dx2iall(ix2)*x%dx2all(ix2+1))+gradSigH3(ix2,ix3)/(x%dx2all(ix2)+x%dx2all(ix2+1))    !static

        coeff=Cmh2(ix2+1,ix3)/(dt*x%dx2iall(ix2)*x%dx2all(ix2+1))
        M(ient)=M(ient)+coeff    !pol. time deriv. terms
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2+1,ix3)    !BC's and pol. time deriv.

        coeff=-1d0*Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+2)*x%dx2iall(ix2+1)) )+ &
              (-1d0)*Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)*x%dx2iall(ix2+1)) )
        M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

        coeff=-1d0*Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )+ &
              (-1d0)*Cm(ix2,ix3-1)*v2(ix2,ix3-1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)+x%dx3all(ix3))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        M(ient)=M(ient)+coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

        ient=ient+1


        !ix2+2,ix3 grid point
        coeff=Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+2)*x%dx2iall(ix2+1)) )
        if (ix2==lx2-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx2(ix3)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi+2

          M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

          ient=ient+1
        end if


        !ix2-2,ix3+1 grid point
        coeff=Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)+x%dx2all(ix2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        if (ix2==2) then
          b(iPhi)=b(iPhi)-coeff*Vminx2(ix3+1)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi+lx2-2

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !!!ix2,ix3+1 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi+lx2

        M(ient)=SigPh3(ix2,ix3+1)/(x%dx3iall(ix3)*x%dx3all(ix3+1))-gradSigH2(ix2,ix3)/(x%dx3all(ix3)+x%dx3all(ix3+1))    !static

        coeff=Cmh3(ix2,ix3+1)/(dt*x%dx3iall(ix3)*x%dx3all(ix3+1))
        M(ient)=M(ient)+coeff    !pol. time deriv.
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2,ix3+1)    !BC's and pol. time deriv.

        coeff=-1d0*Cm(ix2,ix3+1)*v3(ix2,ix3+1)/( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+2)*x%dx3iall(ix3+1)) )+ &
              (-1d0)*Cm(ix2,ix3+1)*v3(ix2,ix3+1)/( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)*x%dx3iall(ix3+1)) )
        M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

        coeff=-1d0*Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)+x%dx2all(ix2+2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )+ &
              (-1d0)*Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)+x%dx2all(ix2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        M(ient)=M(ient)+coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

        ient=ient+1


        !ix2+2,ix3+1 grid point
        coeff=Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)+x%dx2all(ix2+2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        if (ix2==lx2-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx2(ix3+1)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi+lx2+2

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !ix2-1,ix3+2 grid point
        coeff=-1d0*Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
              ( (x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        if (ix3==lx3-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx3(ix2-1)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi+2*lx2-1

          M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

          ient=ient+1
        end if


        !ix2,ix3+2 grid point
        coeff=Cm(ix2,ix3+1)*v3(ix2,ix3+1)/( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+2)*x%dx3iall(ix3+1)) )
        if (ix3==lx3-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx3(ix2)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi+2*lx2

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !ix2+1,ix3+2 grid point
        coeff=Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        if (ix3==lx3-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx3(ix2+1)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi+2*lx2+1

          M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

          ient=ient+1
        end if
      end if
    end do
  end do
end if


!FIRE UP MUMPS
if (myid == 0) then
  print *, 'Filled ',ient-1,' matrix entries.  Initializing MUMPS...'
end if
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!LOAD OUR PROBLEM
if ( myid==0 ) then
  mumps_par%N=lPhi
  mumps_par%NZ=lent
  allocate( mumps_par%IRN ( mumps_par%NZ ) )
  allocate( mumps_par%JCN ( mumps_par%NZ ) )
  allocate( mumps_par%A( mumps_par%NZ ) )
  allocate( mumps_par%RHS ( mumps_par%N  ) )
  mumps_par%IRN=ir
  mumps_par%JCN=ic
  mumps_par%A=M
  mumps_par%RHS=b
  deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

  if (perflag .and. it/=1) then       !used cached permutation
    allocate(mumps_par%PERM_IN(mumps_par%N))
    mumps_par%PERM_IN=mumps_perm
    mumps_par%ICNTL(7)=1
  end if

  !may solve some memory allocation issues, uncomment if MUMPS throws errors
  !about not having enough memory
  !mumps_par%ICNTL(14)=50
end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
if ( myid==0 ) then
  print *, 'Now organizing results...'

  if (perflag .and. it==1) then
    allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
    mumps_perm=mumps_par%SYM_PERM
  end if

  elliptic2D_pol_conv_curv=reshape(mumps_par%RHS,[lx2,lx3])

  print *, 'Now attempting deallocations...'

  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
end if

mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif

end function elliptic2D_pol_conv_curv


function elliptic2D_pol_conv_curv_periodic2(srcterm,SigP,SigH,Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,x,Phi0,perflag,it)

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

real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: SigPh2
real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: SigPh3
!    real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: SigHh2
!    real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: SigHh3
real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: Cmh2
real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: Cmh3
real(wp), dimension(1:size(SigP,1),1:size(SigP,2)+1) :: gradSigH2,gradSigH3

real(wp) :: coeff    !coefficient for calculating polarization terms
integer :: ix2,ix3,lx2,lx3    !this overwrites the values stored in the grid module, which is fine, but perhaps redundant
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
real(wp) :: tstart,tfin

#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#elif REALBITS==64
type (DMUMPS_STRUC) mumps_par
#else
error stop "realbits must be 32 or 64"
#endif

integer :: lcount,ix2tmp,ix3tmp

real(wp), dimension(size(SigP,1),size(SigP,2)) :: tmpresults

real(wp), dimension(size(SigP,1),size(SigP,2)) :: elliptic2D_pol_conv_curv_periodic2


!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
if (myid==0) then
  lx2=size(SigP,1)    !these are full-grid sizes since grid module globals are not in scope
  lx3=size(SigP,2)
  lPhi=lx2*lx3
!      lent=5*(lx2-2)*(lx3-2)+2*lx2+2*(lx3-2)    !static model
!      lent=17*(lx2-2)*(lx3A-2)+2*lx2+2*(lx3-2)-3*2*(lx2-2)-3*2*(lx3-2)    !interior+boundary-x3_adj-x2_adj.  Note that are 3 sets of entries for each adjacent point.  The x3 adjacent points do not need to be removed in teh case of periodic boundary conditions.  This is technicall correct, but a bit misleading, I think.  Shouldn't it be lent=17*(lx2-2)*(lx3-2)+2*(lx2-2)+2*lx3-3*2*(lx2-2)-3*2*(lx3-2)???????
!      lent=17*(lx2-2)*(lx3+1-2)+2*(lx3+1)+3*(lx2-2)-3*2*(lx3+1-2)    !true interior with x3 boundaries which are not treated as interior in periodic solves  + add x2 boundaries (note that these are now size lx3+1) + 3 entries for each x3 ghost cell that we are adding - x2_adj (x2 is not periodici, two sets of three points each, note again the larger x3 size as compared to aperiodic solutions).
  lent=17*(lx2-2)*(lx3-2)+2*(lx3)+17*2*(lx2-2)-3*2*lx3    !true interior with x3 boundaries which are not treated as interior in periodic solves  + add x2 boundaries (note that these are now size lx3+1) + x3 edge cells (treated here as interior) - x2_adj (x2 is not periodici, two sets of three points each, note again the larger x3 size as compared to aperiodic solutions).
  allocate(ir(lent),ic(lent),M(lent),b(lPhi))

  print *, 'MUMPS will attempt a solve of size:  ',lx2,lx3
  print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


  !NOTE THAT THESE NEED TO BE PERIODIC IN X3
  SigPh2(1,:)=0d0
  SigPh2(2:lx2,:)=0.5d0*(SigP(1:lx2-1,:)+SigP(2:lx2,:))
  SigPh3(:,1)=0.5d0*(SigP(:,lx3)+SigP(:,1))    !needs to be left interface value so average of first and last grid point
  SigPh3(:,2:lx3)=0.5d0*(SigP(:,1:lx3-1)+SigP(:,2:lx3))
  Cmh2(1,:)=0d0
  Cmh2(2:lx2,:)=0.5d0*(Cm(1:lx2-1,:)+Cm(2:lx2,:))
  Cmh3(:,1)=0.5d0*(Cm(:,lx3)+Cm(:,1))
  Cmh3(:,2:lx3)=0.5d0*(Cm(:,1:lx3-1)+Cm(:,2:lx3))



  !ZZZ - THESE NEED TO BE CHANGED INTO CIRCULAR/PERIODIC DERIVATIVES FOR THE X3 DIRECTION

  gradSigH2=grad2D1_curv_alt(SigH,x,1,lx2)   !note the alt since we need to use dx2 as differential...  Tricky bug/feature
  gradSigH3=grad2D3_curv_periodic(SigH,x,1,lx3)    !circular difference


  !------------------------------------------------------------
  !-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
  !-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
  !-------GUIDE).
  !------------------------------------------------------------

  print *, 'Loading up matrix entries...'

  !LOAD UP MATRIX ELEMENTS
  lcount=0
  M(:)=0d0
  b=pack(srcterm,.true.)           !boundaries overwritten later, polarization terms also added later.
  ient=1
  do ix3=1,lx3    !note that we have one extra ghost cell now to deal with due to the way we've chosen to implement periodic boundary conditions
    do ix2=1,lx2
      iPhi=lx2*(ix3-1)+ix2     !linear index referencing Phi(ix2,ix3) as a column vector.  Also row of big matrix

      if (ix2==1) then          !BOTTOM GRID POINTS + CORNER
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1d0
        b(iPhi)=Vminx2(ix3)
        ient=ient+1
      elseif (ix2==lx2) then    !TOP GRID POINTS + CORNER
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1d0
        b(iPhi)=Vmaxx2(ix3)
        ient=ient+1
      else                      !TREAT AS AN INTERIOR LOCATION, THIS INCLUDE X3 EDGES NOW SINCE PERIODIC,CIRCULANT


      !ZZZ - NEED TO WRAP INDICES AROUND:  X 1) matrix row/column entries; 2) references to dx3i*(anything but ix3);  3) references to conductances/bcs/etc.

        !ix2-1,ix3-2 grid point
        coeff=-1d0*Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v2(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)+x%dx3all(ix3))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        ir(ient)=iPhi
        ix3tmp=mod(ix3-2-1+lx3,lx3)+1
        ix2tmp=ix2-1
        ic(ient)=lx2*(ix3tmp-1)+ix2tmp
!              ic(ient)=iPhi-2*lx2-1+lPhi-1   !add the grid size to wrap the index around (would be negative otherwise), -1 at end because the last grid opint is actually the same as the first for our implementation
!            else
!              ic(ient)=iPhi-2*lx2-1
!            end if
        M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
        ient=ient+1


        !ix2,ix3-2 grid point
        coeff=-1d0*Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v3(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
                 ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)* &
                 x%dx3iall(mod(ix3-1-1+lx3,lx3)+1)) )
        ir(ient)=iPhi
        ix3tmp=mod(ix3-2-1+lx3,lx3)+1
        ix2tmp=ix2
        ic(ient)=lx2*(ix3tmp-1)+ix2tmp
!              ic(ient)=iPhi-2*lx2+lPhi-1    !again add grid size to wrap to end
!            else
!              ic(ient)=iPhi-2*lx2
!            end if
        M(ient)=coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term
        ient=ient+1


        !ix2+1,ix3-2 grid point
        coeff=Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v2(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)+x%dx3all(ix3))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        ir(ient)=iPhi
        ix3tmp=mod(ix3-2-1+lx3,lx3)+1
        ix2tmp=ix2+1
        ic(ient)=lx2*(ix3tmp-1)+ix2tmp
        M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
        ient=ient+1


        !ix2-2,ix3-1
        coeff=-1d0*Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)+x%dx2all(ix2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        if (ix2==2) then
          b(iPhi)=b(iPhi)-coeff*Vminx2(mod(ix3-1-1+lx3,lx3)+1)
        else
          ir(ient)=iPhi
          ix3tmp=mod(ix3-1-1+lx3,lx3)+1
          ix2tmp=ix2-2
          ic(ient)=lx2*(ix3tmp-1)+ix2tmp
!              ic(ient)=iPhi-lx2-2

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !ix2,ix3-1 grid point in ix2,ix3 equation
        ir(ient)=iPhi
!            ic(ient)=iPhi-lx2
        ix3tmp=mod(ix3-1-1+lx3,lx3)+1
        ix2tmp=ix2
        ic(ient)=lx2*(ix3tmp-1)+ix2tmp

        M(ient)=SigPh3(ix2,ix3)/(x%dx3iall(ix3)*x%dx3all(ix3))+gradSigH2(ix2,ix3)/(x%dx3all(ix3)+x%dx3all(ix3+1))    !static terms

        coeff=Cmh3(ix2,ix3)/(dt*x%dx3iall(ix3)*x%dx3all(ix3))
        M(ient)=M(ient)+coeff   !polarization time derivative terms
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2,mod(ix3-1-1+lx3,lx3)+1)   !add in polarziation terms that include previous time step potential at this grid point

        coeff=Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v3(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
                 ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3)*x%dx3iall(mod(ix3-1-1+lx3,lx3)+1)) )+ &
              Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v3(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
                 ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)*x%dx3iall(mod(ix3-1-1+lx3,lx3)+1)) )
        M(ient)=M(ient)+coeff   !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

        coeff=Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)+x%dx2all(ix2+2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )+ &
              Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)+x%dx2all(ix2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        M(ient)=M(ient)+coeff   !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

        ient=ient+1


        !ix2+2,ix3-1 grid point
        coeff=-1d0*Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
              ( (x%dx2all(ix2+1)+x%dx2all(ix2+2))*(x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        if (ix2==lx2-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx2(mod(ix3-1-1+lx3,lx3)+1)
        else
          ir(ient)=iPhi
!              ic(ient)=iPhi-lx2+2
          ix3tmp=mod(ix3-1-1+lx3,lx3)+1
          ix2tmp=ix2+2
          ic(ient)=lx2*(ix3tmp-1)+ix2tmp

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !ix2-2,ix3 grid point
        coeff=-1d0*Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)*x%dx2iall(ix2-1)) )
        if (ix2==2) then
          b(iPhi)=b(iPhi)-coeff*Vminx2(ix3)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-2

          M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

          ient=ient+1
        end if


        !ix2-1,ix3 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi-1

        M(ient)=SigPh2(ix2,ix3)/(x%dx2iall(ix2)*x%dx2all(ix2))-gradSigH3(ix2,ix3)/(x%dx2all(ix2)+x%dx2all(ix2+1))    !static

        coeff=Cmh2(ix2,ix3)/(dt*x%dx2iall(ix2)*x%dx2all(ix2))
        M(ient)=M(ient)+coeff    !pol. time deriv.
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2-1,ix3)    !BC's and pol. time deriv.

        coeff=Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2)*x%dx2iall(ix2-1)) )+ &
              Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)*x%dx2iall(ix2-1)) )
        M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

        coeff=Cm(ix2,mod(ix3+1-1,lx3)+1)*v2(ix2,mod(ix3+1-1,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )+ &
              Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v2(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)+x%dx3all(ix3))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        M(ient)=M(ient)+coeff     !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

        ient=ient+1


        !ix2,ix3 grid point (main diagonal)
        ir(ient)=iPhi
        ic(ient)=iPhi

        M(ient)=-1d0*SigPh2(ix2+1,ix3)/(x%dx2iall(ix2)*x%dx2all(ix2+1)) &
                -1d0*SigPh2(ix2,ix3)/(x%dx2iall(ix2)*x%dx2all(ix2)) &
                -1d0*SigPh3(ix2,mod(ix3+1-1,lx3)+1)/(x%dx3iall(ix3)*x%dx3all(ix3+1)) &
                -1d0*SigPh3(ix2,ix3)/(x%dx3iall(ix3)*x%dx3all(ix3))    !static

        coeff=-1d0*Cmh2(ix2+1,ix3)/(dt*x%dx2iall(ix2)*x%dx2all(ix2+1)) &
              -1d0*Cmh2(ix2,ix3)/(dt*x%dx2iall(ix2)*x%dx2all(ix2)) &
              -1d0*Cmh3(ix2,mod(ix3+1-1,lx3)+1)/(dt*x%dx3iall(ix3)*x%dx3all(ix3+1)) &
              -1d0*Cmh3(ix2,ix3)/(dt*x%dx3iall(ix3)*x%dx3all(ix3))
        M(ient)=M(ient)+coeff    !pol. time deriv.
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2,ix3)    !BC's and pol. time deriv.

        coeff=Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)*x%dx2iall(ix2+1)) )+ &
              (-1d0)*Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2)*x%dx2iall(ix2-1)) )
        M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

        coeff=Cm(ix2,mod(ix3+1-1,lx3)+1)*v3(ix2,mod(ix3+1-1,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)*x%dx3all(ix3+1)) )+ &
              (-1d0)*Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v3(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3)*x%dx3iall(mod(ix3-1-1+lx3,lx3)+1)) )
        M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

        ient=ient+1


        !ix2+1,ix3 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi+1

        M(ient)=SigPh2(ix2+1,ix3)/(x%dx2iall(ix2)*x%dx2all(ix2+1))+gradSigH3(ix2,ix3)/(x%dx2all(ix2)+x%dx2all(ix2+1))    !static

        coeff=Cmh2(ix2+1,ix3)/(dt*x%dx2iall(ix2)*x%dx2all(ix2+1))
        M(ient)=M(ient)+coeff    !pol. time deriv. terms
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2+1,ix3)    !BC's and pol. time deriv.

        coeff=-1d0*Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+2)*x%dx2iall(ix2+1)) )+ &
              (-1d0)*Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)*x%dx2iall(ix2+1)) )
        M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

        coeff=-1d0*Cm(ix2,mod(ix3+1-1,lx3)+1)*v2(ix2,mod(ix3+1-1,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )+ &
              (-1d0)*Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v2(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3-1)+x%dx3all(ix3))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        M(ient)=M(ient)+coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

        ient=ient+1


        !ix2+2,ix3 grid point
        coeff=Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+2)*x%dx2iall(ix2+1)) )
        if (ix2==lx2-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx2(ix3)
        else
          ir(ient)=iPhi
          ic(ient)=iPhi+2

          M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

          ient=ient+1
        end if


        !ix2-2,ix3+1 grid point
        coeff=Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)+x%dx2all(ix2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        if (ix2==2) then
          b(iPhi)=b(iPhi)-coeff*Vminx2(mod(ix3+1-1,lx3)+1)
        else
          ir(ient)=iPhi
!              ic(ient)=iPhi+lx2-2
          ix3tmp=mod(ix3+1-1,lx3)+1
          ix2tmp=ix2-2
          ic(ient)=lx2*(ix3tmp-1)+ix2tmp

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !ix2,ix3+1 grid point
        ir(ient)=iPhi
!            ic(ient)=iPhi+lx2
        ix3tmp=mod(ix3+1-1,lx3)+1
        ix2tmp=ix2
        ic(ient)=lx2*(ix3tmp-1)+ix2tmp

        M(ient)=SigPh3(ix2,mod(ix3+1-1,lx3)+1)/ &
                (x%dx3iall(ix3)*x%dx3all(ix3+1))-gradSigH2(ix2,ix3)/(x%dx3all(ix3)+x%dx3all(ix3+1))    !static

        coeff=Cmh3(ix2,mod(ix3+1-1,lx3)+1)/(dt*x%dx3iall(ix3)*x%dx3all(ix3+1))
        M(ient)=M(ient)+coeff    !pol. time deriv.
        b(iPhi)=b(iPhi)+coeff*Phi0(ix2,mod(ix3+1-1,lx3)+1)    !BC's and pol. time deriv.

        coeff=-1d0*Cm(ix2,mod(ix3+1-1,lx3)+1)*v3(ix2,mod(ix3+1-1,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+2)*x%dx3iall(mod(ix3+1-1,lx3)+1)) )+ &
              (-1d0)*Cm(ix2,mod(ix3+1-1,lx3)+1)*v3(ix2,mod(ix3+1-1,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)*x%dx3iall(mod(ix3+1-1,lx3)+1)) )
        M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

        coeff=-1d0*Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)+x%dx2all(ix2+2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )+ &
              (-1d0)*Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2-1)+x%dx2all(ix2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        M(ient)=M(ient)+coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

        ient=ient+1


        !ix2+2,ix3+1 grid point
        coeff=Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
              ( (x%dx2all(ix2)+x%dx2all(ix2+1))*(x%dx2all(ix2+1)+x%dx2all(ix2+2))*(x%dx3all(ix3)+x%dx3all(ix3+1)) )
        if (ix2==lx2-1) then
          b(iPhi)=b(iPhi)-coeff*Vmaxx2(mod(ix3+1-1,lx3)+1)
        else
          ir(ient)=iPhi
!              ic(ient)=iPhi+lx2+2
          ix3tmp=mod(ix3+1-1,lx3)+1
          ix2tmp=ix2+2
          ic(ient)=lx2*(ix3tmp-1)+ix2tmp

          M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

          ient=ient+1
        end if


        !ix2-1,ix3+2 grid point
!            coeff=-1d0*Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
!                  ( (x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        coeff=-1d0*Cm(ix2,mod(ix3+1-1,lx3)+1)*v2(ix2,mod(ix3+1-1,lx3)+1)/ &      !mods to wrap indices around
              ( (x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        ir(ient)=iPhi
!            if (ix3>=lx3-1) then    !this needs to also handle the case where ix3=lx3!!!  Likewise for statements that follow...
          ix3tmp=mod(ix3+2-1,lx3)+1
          ix2tmp=ix2-1
          ic(ient)=lx2*(ix3tmp-1)+ix2tmp
!              ic(ient)=iPhi+2*lx2-1-lPhi+1    !wrap to beginning
!            else
!              ic(ient)=iPhi+2*lx2-1
!            end if
        M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
        ient=ient+1


        !ix2,ix3+2 grid point
!            coeff=Cm(ix2,ix3+1)*v3(ix2,ix3+1)/ &
!                  ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+2)*x%dx3iall(ix3+1)) )
        coeff=Cm(ix2,mod(ix3+1-1,lx3)+1)*v3(ix2,mod(ix3+1-1,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+2)*x%dx3iall(mod(ix3+1-1,lx3)+1)) )
        ir(ient)=iPhi
!            if (ix3>=lx3-1) then
          ix3tmp=mod(ix3+2-1,lx3)+1
          ix2tmp=ix2
          ic(ient)=lx2*(ix3tmp-1)+ix2tmp
!              ic(ient)=iPhi+2*lx2-lPhi+1    !substract grid size to wrap around to the beginning
!            else
!              ic(ient)=iPhi+2*lx2
!            end if
        M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )
        ient=ient+1


        !ix2+1,ix3+2 grid point
!            coeff=Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
!                  ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        coeff=Cm(ix2,mod(ix3+1-1,lx3)+1)*v2(ix2,mod(ix3+1-1,lx3)+1)/ &
              ( (x%dx3all(ix3)+x%dx3all(ix3+1))*(x%dx3all(ix3+1)+x%dx3all(ix3+2))*(x%dx2all(ix2)+x%dx2all(ix2+1)) )
        ir(ient)=iPhi
!            if (ix3>=lx3-1) then     !this should actually work for any case...
          ix3tmp=mod(ix3+2-1,lx3)+1
          ix2tmp=ix2+1
          ic(ient)=lx2*(ix3tmp-1)+ix2tmp
!              ic(ient)=iPhi+2*lx2+1-lPhi+1
!            else
!              ic(ient)=iPhi+2*lx2+1
!            end if
        M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
        ient=ient+1
      end if
    end do
  end do
end if


!FIRE UP MUMPS
if (myid == 0) then
  print *,  'Debug count:  ',lcount
  print *, 'Filled ',ient-1,' out of ',lent,' matrix entries for solving ',iPhi,' of ',lPhi, &
               ' unknowns.  Initializing MUMPS...'
  if (ient-1 /= lent) then
    error stop 'Incorrect number of matrix entries filled in potential solve!!!'
  end if
end if
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!LOAD OUR PROBLEM
if ( myid==0 ) then
  mumps_par%N=lPhi
  mumps_par%NZ=lent
  allocate( mumps_par%IRN ( mumps_par%NZ ) )
  allocate( mumps_par%JCN ( mumps_par%NZ ) )
  allocate( mumps_par%A( mumps_par%NZ ) )
  allocate( mumps_par%RHS ( mumps_par%N  ) )
  mumps_par%IRN=ir
  mumps_par%JCN=ic
  mumps_par%A=M
  mumps_par%RHS=b
  deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

  if (perflag .and. it/=1) then       !used cached permutation
    allocate(mumps_par%PERM_IN(mumps_par%N))
    mumps_par%PERM_IN=mumps_perm
    mumps_par%ICNTL(7)=1
  end if

  !mumps_par%ICNTL(14)=50
end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
if ( myid==0 ) then
  print *, 'Now organizing results...'

  if (perflag .and. it==1) then
    allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
    mumps_perm=mumps_par%SYM_PERM
  end if

 !IF WE HAVE DONE A PERIODIC SOLVE, THE LAST GRID POINT NEEDS TO BE IGNORED WHEN WE RESHAPE THE POTENTIAL ARRAY.

  tmpresults=reshape(mumps_par%RHS,[lx2,lx3])
  elliptic2D_pol_conv_curv_periodic2=tmpresults(1:lx2,1:lx3)    !sort of superfluous now that hte solve size is the same as the grid

  print *, 'Now attempting deallocations...'

  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
end if

mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif

end function elliptic2D_pol_conv_curv_periodic2


function elliptic2D_nonint_curv(srcterm,sig0,sigP,Vminx1,Vmaxx1,Vminx3,Vmaxx3,x,flagdirich,perflag,it)

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

real(wp), dimension(1:size(sig0,1),1:size(sig0,3)) :: sig0h1
real(wp), dimension(1:size(sigP,1),1:size(sigP,3)) :: sigPh3

integer :: ix1,ix3,lx1,lx3
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
real(wp) :: tstart,tfin

#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#elif REALBITS==64
type (DMUMPS_STRUC) mumps_par
#else
error stop "realbits must be 32 or 64"
#endif

real(wp), dimension(size(sig0,1),1,size(sig0,3)) :: elliptic2D_nonint_curv


!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
if (myid==0) then
  lx1=size(sig0,1)
  lx3=size(sig0,3)
  lPhi=lx1*lx3
  if (flagdirich==0) then
    lent=5*(lx1-2)*(lx3-2)+2*lx1+2*(lx3-2)+2*lx3    !first +1 for Neumann bottom, second for Neumann top
  else
!        lent=5*(lx1-2)*(lx3-2)+2*lx1+2*(lx3-2)+1    !first +1 for Neumann bottom
    lent=5*(lx1-2)*(lx3-2)+2*(lx1-2)+2*lx3+lx3
  end if
  allocate(ir(lent),ic(lent),M(lent),b(lPhi))

  print *, 'MUMPS will attempt a solve of size:  ',lx1,lx3
  print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


  !AVERAGE CONDUCTANCES TO THE GRID INTERFACE POINTS
  sig0h1(1,:)=0d0
  sig0h1(2:lx1,:)=0.5d0*(sig0(1:lx1-1,1,:)+sig0(2:lx1,1,:))
  sigPh3(:,1)=0d0
  sigPh3(:,2:lx3)=0.5d0*(sigP(:,1,1:lx3-1)+sigP(:,1,2:lx3))

  !------------------------------------------------------------
  !-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
  !-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
  !-------GUIDE).
  !------------------------------------------------------------

  !LOAD UP MATRIX ELEMENTS
  M(:)=0d0
  b=pack(srcterm,.true.)           !boundaries overwritten later
  ient=1
  do ix3=1,lx3
    do ix1=1,lx1
      iPhi=lx1*(ix3-1)+ix1     !linear index referencing Phi(ix1,ix3) as a column vector.  Also row of big matrix

      if (ix1==1) then          !(LOGICAL) BOTTOM GRID POINTS + CORNER, USE NEUMANN HERE, PRESUMABLY ZERO.  ZZZ - potential problem here if we have inverted grid...
        if (gridflag/=1) then
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=-1d0
          ient=ient+1
          ir(ient)=iPhi
          ic(ient)=iPhi+1
          M(ient)=1d0
!            b(iPhi)=Vminx1(ix3)
          b(iPhi)=0d0    !force bottom current to zero
          ient=ient+1
        else
          if (flagdirich/=0) then    !ZZZ - need to check non-inverted???
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1d0
            b(iPhi)=Vmaxx1(1,ix3)
            ient=ient+1
          else
            ir(ient)=iPhi
            ic(ient)=iPhi-1
            M(ient)=-1d0/x%dx1(lx1)
            ient=ient+1
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1d0/x%dx1(lx1)
            b(iPhi)=Vmaxx1(1,ix3)
            ient=ient+1
          end if
        end if
      elseif (ix1==lx1) then    !(LOGICAL) TOP GRID POINTS + CORNER
        if (gridflag/=1) then    !non-inverted?
          if (flagdirich/=0) then    !ZZZ - need to check non-inverted???
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1d0
            b(iPhi)=Vmaxx1(1,ix3)
            ient=ient+1
          else
            ir(ient)=iPhi
            ic(ient)=iPhi-1
            M(ient)=-1d0/x%dx1(lx1)
            ient=ient+1
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1d0/x%dx1(lx1)
            b(iPhi)=Vmaxx1(1,ix3)
            ient=ient+1
          end if
        else
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=-1d0
          ient=ient+1
          ir(ient)=iPhi
          ic(ient)=iPhi+1
          M(ient)=1d0
!            b(iPhi)=Vminx1(ix3)
          b(iPhi)=0d0    !force bottom current to zero
          ient=ient+1
        end if
      elseif (ix3==1) then      !LEFT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vminx3(ix1,1)
        ient=ient+1
      elseif (ix3==lx3) then    !RIGHT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vmaxx3(ix1,1)
        ient=ient+1
      else                      !INTERIOR
        !ix1,ix3-1 grid point in ix1,ix3 equation
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1
        M(ient)=sigPh3(ix1,ix3)/(x%dx3iall(ix3)*x%dx3all(ix3))
        ient=ient+1

        !ix1-1,ix3 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi-1
        M(ient)=sig0h1(ix1,ix3)/(x%dx1i(ix1)*x%dx1(ix1))
        ient=ient+1

        !ix1,ix3 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-1d0*sig0h1(ix1+1,ix3)/(x%dx1i(ix1)*x%dx1(ix1+1)) &
                -1d0*sig0h1(ix1,ix3)/(x%dx1i(ix1)*x%dx1(ix1)) &
                -1d0*sigPh3(ix1,ix3+1)/(x%dx3iall(ix3)*x%dx3all(ix3+1)) &
                -1d0*sigPh3(ix1,ix3)/(x%dx3iall(ix3)*x%dx3all(ix3))
        ient=ient+1

        !ix1+1,ix3 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient)=sig0h1(ix1+1,ix3)/(x%dx1i(ix1)*x%dx1(ix1+1))
        ient=ient+1

        !ix1,ix3+1 grid point
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1
        M(ient)=sigPh3(ix1,ix3+1)/(x%dx3iall(ix3)*x%dx3all(ix3+1))
        ient=ient+1
      end if
    end do
  end do
end if
print *, 'Number of entries used:  ',ient-1


!FIRE UP MUMPS
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!LOAD OUR PROBLEM
if ( myid==0 ) then
  mumps_par%N=lPhi
  mumps_par%NZ=lent
  allocate( mumps_par%IRN ( mumps_par%NZ ) )
  allocate( mumps_par%JCN ( mumps_par%NZ ) )
  allocate( mumps_par%A( mumps_par%NZ ) )
  allocate( mumps_par%RHS ( mumps_par%N  ) )
  mumps_par%IRN=ir
  mumps_par%JCN=ic
  mumps_par%A=M
  mumps_par%RHS=b
  deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

  if (perflag .and. it/=1) then       !used cached permutation
    print *, 'Using a previously stored permutation'
    allocate(mumps_par%PERM_IN(mumps_par%N))
    mumps_par%PERM_IN=mumps_perm
    mumps_par%ICNTL(7)=1
  end if

  !may solve some memory allocation issues, uncomment if MUMPS throws errors
  !about not having enough memory
  !mumps_par%ICNTL(14)=50
end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors).  WOW - this halves execution
! time on some big 2048*2048 solves!!!
if ( myid==0 ) then
  print *, 'Now organizing results...'

  if (perflag .and. it==1) then
    print *, 'Storing ordering for future time step use...'
    allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
    mumps_perm=mumps_par%SYM_PERM
  end if

  elliptic2D_nonint_curv=reshape(mumps_par%RHS,[lx1,1,lx3])

  print *, 'Now attempting deallocations...'

  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
end if

mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif

end function elliptic2D_nonint_curv


function poisson2D(rho,Vminx1,Vmaxx1,Vminx2,Vmaxx2,dx1,perflag)

!------------------------------------------------------------
!-------SOLVE POISSONS'S EQUATION IN 2D USING MUMPS.  THIS MAY
!-------BE USED BY SOME OF THE UNIT TEST PROGRAMS, BUT I FORGET...
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: rho
real(wp), dimension(:), intent(in) :: Vminx1,Vmaxx1
real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
real(wp), intent(in) :: dx1
!    real(wp), dimension(:), intent(in) :: dx2
logical, intent(in) :: perflag

integer :: ix1,ix2,lx1,lx2
integer :: lPhi, lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
real(wp) :: tstart,tfin

#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#elif REALBITS==64
type (DMUMPS_STRUC) mumps_par
#else
error stop "realbits must be 32 or 64"
#endif
real(wp), dimension(size(rho,1),size(rho,2)) :: poisson2D


if (myid==0) then
!------------------------------------------------------------
!-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
!-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
!-------GUIDE).
!------------------------------------------------------------
lx1=size(rho,1)
lx2=size(rho,2)
lPhi=lx1*lx2
lent=5*(lx1-2)*(lx2-2)+2*lx1+2*(lx2-2)
allocate(ir(lent),ic(lent),M(lent),b(lPhi))


!LOAD UP MATRIX ELEMENTS
M(:)=0.0
b=pack(rho,.true.)           !boundaries overwritten later
ient=1
do ix2=1,lx2
  do ix1=1,lx1
    iPhi=lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix2) as a column vector.  Also row of big matrix

    if (ix1==1) then          !BOTTOM GRID POINTS + CORNER
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1.0
      b(iPhi)=Vminx1(ix2)
      ient=ient+1
    elseif (ix1==lx1) then    !TOP GRID POINTS + CORNER
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1.0
      b(iPhi)=Vmaxx1(ix2)
      ient=ient+1
    elseif (ix2==1) then      !LEFT BOUNDARY
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1.0
      b(iPhi)=Vminx2(ix1)
      ient=ient+1
    elseif (ix2==lx2) then    !RIGHT BOUNDARY
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1.0
      b(iPhi)=Vmaxx2(ix1)
      ient=ient+1
    else                      !INTERIOR
      !ix1,ix2-1 grid point in ix1,ix2 equation
      ir(ient)=iPhi
      ic(ient)=iPhi-lx1
      M(ient)=1.0
      ient=ient+1

      !ix1-1,ix2 grid point
      ir(ient)=iPhi
      ic(ient)=iPhi-1
      M(ient)=1.0
      ient=ient+1

      !ix1,ix2 grid point
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=-4.0
      ient=ient+1

      !ix1+1,ix2 grid point
      ir(ient)=iPhi
      ic(ient)=iPhi+1
      M(ient)=1.0
      ient=ient+1

      !ix1,ix2+1 grid point
      ir(ient)=iPhi
      ic(ient)=iPhi+lx1
      M(ient)=1.0
      ient=ient+1
    end if
  end do
end do


!CORRECT FOR DX /= 1
b=b*dx1**2
end if

!FIRE UP MUMPS
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!LOAD OUR PROBLEM (ROOT ONLY)
if ( myid == 0 ) then
  mumps_par%N=lPhi
  mumps_par%NZ=lent
  allocate( mumps_par%IRN ( mumps_par%NZ ) )
  allocate( mumps_par%JCN ( mumps_par%NZ ) )
  allocate( mumps_par%A( mumps_par%NZ ) )
  allocate( mumps_par%RHS ( mumps_par%N  ) )
  mumps_par%IRN=ir
  mumps_par%JCN=ic
  mumps_par%A=M
  mumps_par%RHS=b
  deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

  if (perflag) then       !used cached permutation
    allocate(mumps_par%PERM_IN(mumps_par%N))
    mumps_par%PERM_IN=mumps_perm
    mumps_par%ICNTL(7)=1
  end if
end if


!SOLVE OUR PROBLEM (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
call cpu_time(tstart)
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif
call cpu_time(tfin)
print *, 'Solve took ',tfin-tstart,' seconds...'


!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
if ( myid == 0 ) then
  mumps_perm=mumps_par%SYM_PERM
  poisson2D=reshape(mumps_par%RHS/dx1**2,[lx1,lx2])

  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
end if

mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif

end function poisson2D


subroutine elliptic_workers()

!------------------------------------------------------------
!-------A FN. THAT ALLOWS WORKERS TO ENTER MUMPS SOLVES
!------------------------------------------------------------
#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#elif REALBITS==64
type (DMUMPS_STRUC) mumps_par
#else
error stop "realbits must be 32 or 64"
#endif


!FIRE UP MUMPS
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!ROOT WILL LOAD OUR PROBLEM


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#else
error stop "realbits must be 32 or 64"
#endif


!DEALLOCATE STRUCTURES USED BY WORKERS DURING SOLVE
mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#elif REALBITS==64
call DMUMPS(mumps_par)
#endif

end subroutine elliptic_workers

end module potential_mumps
