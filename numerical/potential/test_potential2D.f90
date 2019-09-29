program test_potential2D

!-------------------------------------------------------------------------------
!-------SOLVE LAPLACE'S EQUATION IN 2D USING PDEelliptic, mumps-based libraries
!-------------------------------------------------------------------------------

use mpi
use phys_consts, only: wp,debug,pi
use PDEelliptic, only: elliptic2D_polarization,elliptic2D_cart,elliptic_workers
implicit none

integer, parameter :: lx1=256,lx2=256,lx3=256
integer :: ix1,ix2,ix3

real(wp), dimension(-1:lx1+2) :: x1
real(wp), dimension(-1:lx2+2) :: x2
real(wp), dimension(-1:lx3+2) :: x3
real(wp), dimension(1:lx1+1) :: x1i
real(wp), dimension(1:lx2+2) :: x2i
real(wp), dimension(1:lx3+2) :: x3i
real(wp), dimension(0:lx1+2) :: dx1
real(wp), dimension(0:lx2+2) :: dx2
real(wp), dimension(0:lx3+2) :: dx3
real(wp), dimension(1:lx1) :: dx1i
real(wp), dimension(1:lx2) :: dx2i
real(wp), dimension(1:lx3) :: dx3i

real(wp), dimension(lx3) :: Vminx2,Vmaxx2
real(wp), dimension(lx2) :: Vminx3,Vmaxx3
real(wp), dimension(1,lx3) :: Vminx22,Vmaxx22
real(wp), dimension(lx2,1) :: Vminx32,Vmaxx32
real(wp) :: tstart,tfin
integer :: u,ierr,myid,lid

real(wp), dimension(lx2,lx3) :: Phi,Phi2squeeze,Phitrue,errorMUMPS,errorMUMPS2
real(wp), dimension(lx2,1,lx3) :: Phi2    !FA solver requires different shaped arrays...

real(wp), dimension(lx2,lx3) :: Phi0=1.0_wp,v2=1.0_wp,v3=1.0_wp     !shouldn't be used if D=0
real(wp), dimension(lx2,lx3) :: A=1.0_wp,Ap=1.0_wp,App=0.0_wp,B=0.0_wp,C=0.0_wp,D=0.0_wp
real(wp), dimension(lx2,1,lx3) :: A2=1.0_wp,Ap2=1.0_wp

real(wp), dimension(lx2,lx3) :: srcterm=0.0_wp
real(wp), dimension(lx2,1,lx3) :: srcterm2=0.0_wp
logical :: perflag=.false.     !shouldn't be used
integer :: it=1                !not used
real(wp) :: dt=1.0_wp          !not used
integer :: gridflag=1
integer :: flagdirich=1        !denoting non-inverted grid...


character(*), parameter :: outfile='test_potential2D.dat'


!! mpi starting
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
call mpi_comm_size(MPI_COMM_WORLD,lid,ierr)


!! Set things up to give debug output
debug=.true.


!! Set up grid and compute differences needed for solution of PDE
x1=[ (real(ix1-1,wp)/real(lx1-1,wp), ix1=-1,lx1+2) ]
dx1=x1(0:lx1+2)-x1(-1:lx1+1)
x1i(1:lx1+1)=0.5*(x1(0:lx1)+x1(1:lx1+1))
dx1i=x1i(2:lx1+1)-x1i(1:lx1)

x2=[ (real(ix2-1,wp)/real(lx2-1,wp), ix2=-1,lx2+2) ]
dx2=x2(0:lx2+2)-x2(-1:lx2+1)
x2i(1:lx2+1)=0.5*(x2(0:lx2)+x2(1:lx2+1))
dx2i=x2i(2:lx2+1)-x2i(1:lx2)

x3=[ (real(ix3-1,wp)/real(lx3-1,wp), ix3=-1,lx3+2) ]
dx3=x3(0:lx3+2)-x3(-1:lx3+1)
x3i(1:lx3+1)=0.5*(x3(0:lx3)+x3(1:lx3+1))
dx3i=x3i(2:lx3+1)-x3i(1:lx3)


!! Define boundary conditions for this problem
Vminx2(1:lx3)=0.0_wp
Vmaxx2(1:lx3)=0.0_wp
Vminx3(1:lx2)=0.0_wp
Vmaxx3(1:lx2)=sin(2*pi*x2(1:lx2))

Vminx22(1,1:lx3)=0.0_wp
Vmaxx22(1,1:lx3)=0.0_wp
Vminx32(1:lx2,1)=0.0_wp
Vmaxx32(1:lx2,1)=sin(2*pi*x2(1:lx2))


!! Make the call to PDE elliptic solver library, note the separate calls for root vs. workers
if (myid==0) then
  print*, 'Starting MUMPS solve...'
  Phi=elliptic2D_polarization(srcterm,A,Ap,App,B,C,D,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,dx1, &
                                 dx1i,dx2,dx2i,dx3,dx3i,Phi0,perflag,it)
  Phi2=elliptic2D_cart(srcterm2,A2,Ap2,Vminx22,Vmaxx22,Vminx32,Vmaxx32,dx2,dx2i,dx3,dx3i,flagdirich,perflag,gridflag,it)
  Phi2squeeze(:,:)=Phi2(:,1,:)
  print*, 'MUMPS solve is complete...'
else
  call elliptic_workers()
  call elliptic_workers()    !twice for two calls by root for two different numerical routines
end if


!! Also have root compute an analytical solution and test against numerical
if (myid==0) then
  do ix3=1,lx3
    do ix2=1,lx2
      Phitrue(ix2,ix3)=sinh(2*pi*x3(ix3))/sinh(2*pi)*sin(2*pi*x2(ix2))    !analytical solution for this prolbem (see documentation)
      errorMUMPS(ix2,ix3)=Phi(ix2,ix3)-Phitrue(ix2,ix3)
      errorMUMPS2(ix2,ix3)=Phi2squeeze(ix2,ix3)-Phitrue(ix2,ix3)
    end do
  end do
end if


!! Write some output for visualizations
if (myid==0) then
  print*, 'Numerical solution range:  ',minval(Phi),maxval(Phi)
  print*, 'Analytical solution range:  ',minval(Phitrue),maxval(Phitrue)

  print *,'Root process is writing ',outfile
  open(newunit=u,file=outfile,status='replace')
  write(u,*) lx2
  call writearray(u,x2(1:lx2))
  write(u,*) lx3
  call writearray(u,x3(1:lx3))
  call write2Darray(u,Phi)
  call write2Darray(u,Phi2squeeze)
  call write2Darray(u,Phitrue)
  close(u)

  print*, '1:  Max error over grid:  ',maxval(abs(errorMUMPS))
  print*, '2:  Max error over grid:  ',maxval(abs(errorMUMPS2))
  if (maxval(abs(errorMUMPS))>0.05_wp) error stop '1:  Numerical error too large; check setup/output!!!'
  if (maxval(abs(errorMUMPS2))>0.05_wp) error stop '2:  Numerical error too large; check setup/output!!!'
end if

call mpi_finalize(ierr)


contains

  subroutine writearray(fileunit,array)
    integer, intent(in) :: fileunit
    real(wp), dimension(:), intent(in) :: array

    integer :: k

    do k=1,size(array)
      write(fileunit,*) array(k)
    end do
  end subroutine writearray


  subroutine write2Darray(fileunit,array)
    integer, intent(in) :: fileunit
    real(wp), dimension(:,:), intent(in) :: array

    integer :: k1,k2

    do k1=1,size(array,1)
      write(fileunit,'(f12.6)') (array(k1,k2), k2=1,size(array,2))
    end do
  end subroutine write2Darray

end program test_potential2D
