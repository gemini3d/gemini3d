program test_potential2d
!! SOLVE LAPLACE'S EQUATION IN 2D USING PDEelliptic, mumps-based libraries

use mpi_f08, only : mpi_init, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_finalize
use phys_consts, only: wp,debug,pi
use PDEelliptic, only: elliptic2D_polarization,elliptic2D_cart,elliptic_workers
use h5fortran, only: hdf5_file

implicit none (type, external)

type(hdf5_file) :: hout

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

integer :: ierr,myid,lid

real(wp), allocatable, dimension(:,:) :: Phi,Phi2squeeze,Phitrue,errorMUMPS,errorMUMPS2
real(wp), allocatable, dimension(:,:,:) :: Phi2
!! FA solver requires different shaped arrays...

real(wp), allocatable, dimension(:,:) :: Phi0, v2, v3
!! shouldn't be used if D=0
real(wp), allocatable, dimension(:,:) :: A, Ap, App, B, C, D
real(wp), allocatable, dimension(:,:,:) :: A2, Ap2

real(wp), allocatable, dimension(:,:) :: srcterm
real(wp), allocatable, dimension(:,:,:) :: srcterm2
logical :: perflag=.false.     !shouldn't be used
integer :: it=1                !not used
real(wp) :: dt=1          !not used
integer :: gridflag=1
integer, dimension(4), parameter :: flagsdirich=[1,1,1,1]        !denoting all Dirichlet conditions for test problem

character(4096) :: argv

! --- avoid stack issues by using allocatable()

allocate(Phi(lx2,lx3), Phi2squeeze(lx2,lx3), Phitrue(lx2,lx3), errorMUMPS(lx2,lx3), errorMUMPS2(lx2,lx3))
allocate(Phi2(lx2,1,lx3))

allocate(Phi0(lx2,lx3), v2(lx2,lx3), v3(lx2,lx3))
allocate(A(lx2,lx3), Ap(lx2,lx3), App(lx2,lx3), B(lx2,lx3), C(lx2,lx3), D(lx2,lx3))
Phi0=1
v2=1
v3=1
A=1
Ap=1
App=0
B=0
C=0
D=0

allocate(A2(lx2,1,lx3), Ap2(lx2,1,lx3))
A2=1
Ap2=1

allocate(srcterm(lx2,lx3), srcterm2(lx2,1,lx3))
srcterm=0
srcterm2=0

call mpi_init()
call mpi_comm_rank(MPI_COMM_WORLD,myid)
call mpi_comm_size(MPI_COMM_WORLD,lid)


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
Vminx2(1:lx3)=0
Vmaxx2(1:lx3)=0
Vminx3(1:lx2)=0
Vmaxx3(1:lx2)=sin(2*pi*x2(1:lx2))

Vminx22(1,1:lx3)=0
Vmaxx22(1,1:lx3)=0
Vminx32(1:lx2,1)=0
Vmaxx32(1:lx2,1)=sin(2*pi*x2(1:lx2))


!! Make the call to PDE elliptic solver library, note the separate calls for root vs. workers
if (myid==0) then
  print*, 'Starting MUMPS solve...'
    Phi=elliptic2D_polarization(srcterm,A,Ap,App,B,C,D,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,dx1, &
                                 dx1i,dx2,dx2i,dx3,dx3i,Phi0,perflag,it)
    Phi2=elliptic2D_cart(srcterm2,A2,Ap2,Vminx22,Vmaxx22,Vminx32,Vmaxx32,dx2,dx2i,dx3,dx3i,flagsdirich,perflag,gridflag,it)
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


! Write some output for visualizations
if (myid==0) then
  call get_command_argument(1, argv, status=ierr)
  if(ierr /= 0) error stop 'please specify filename'

  print*, 'Numerical solution range:  ',minval(Phi),maxval(Phi)
  print*, 'Analytical solution range:  ',minval(Phitrue),maxval(Phitrue)

  call hout%open(trim(argv), action="w")
  call hout%write("/lx1", lx1)
  call hout%write("/lx2", lx2)
  call hout%write("/lx3", lx3)
  call hout%write("/x1", x1(1:lx1))
  call hout%write("/x2", x2(1:lx2))
  call hout%write("/x3", x3(1:lx3))
  call hout%write("/Phi", Phi)
  call hout%write("/Phi2squeeze", Phi2squeeze)
  call hout%write("/Phitrue", Phitrue)
  call hout%close()

  print*, '1:  Max error over grid:  ',maxval(abs(errorMUMPS))
  print*, '2:  Max error over grid:  ',maxval(abs(errorMUMPS2))
  if (maxval(abs(errorMUMPS))>0.05_wp) error stop '1:  Numerical error too large; check setup/output!!!'
  if (maxval(abs(errorMUMPS2))>0.05_wp) error stop '2:  Numerical error too large; check setup/output!!!'
end if

call mpi_finalize()

end program

! block
!   integer :: u
!   open(newunit=u, file=outfile, status='replace', action='write')
!   write(u,*) lx2
!   call writearray(u,x2(1:lx2))
!   write(u,*) lx3
!   call writearray(u,x3(1:lx3))
!   call write2Darray(u,Phi)
!   call write2Darray(u,Phi2squeeze)
!   call write2Darray(u,Phitrue)
!   close(u)
! end block

! contains

!   subroutine writearray(fileunit,array)
!     integer, intent(in) :: fileunit
!     real(wp), dimension(:), intent(in) :: array

!     integer :: k

!     do k=1,size(array)
!       write(fileunit,*) array(k)
!     end do
!   end subroutine writearray


!   subroutine write2Darray(fileunit,array)
!     integer, intent(in) :: fileunit
!     real(wp), dimension(:,:), intent(in) :: array

!     integer :: k1,k2

!     do k1=1,size(array,1)
!       write(fileunit,'(f12.6)') (array(k1,k2), k2=1,size(array,2))
!     end do
!   end subroutine write2Darray

! end program test_potential2D
