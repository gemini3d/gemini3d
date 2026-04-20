!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Solve an aurora-like potential problem using MUMPs and GEMINI interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_potential2d_auroral

use mpi_f08, only : mpi_init, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_finalize
use phys_consts, only: wp,debug,pi
use PDEelliptic, only: elliptic2D_static,elliptic_workers
use h5fortran, only: hdf5_file

implicit none (type, external)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type(hdf5_file) :: hout

!! system size
integer, parameter :: lx1=256,lx2=256,lx3=256
integer :: ix1,ix2,ix3
real(wp), parameter :: x2dist=400e3, x3dist=1000e3

!! coordinates
real(wp), dimension(-1:lx1+2) :: x1
real(wp), dimension(-1:lx2+2) :: x2
real(wp), dimension(-1:lx3+2) :: x3
real(wp), dimension(1:lx1+1) :: x1i
real(wp), dimension(1:lx2+1) :: x2i
real(wp), dimension(1:lx3+1) :: x3i
real(wp), dimension(0:lx1+2) :: dx1
real(wp), dimension(0:lx2+2) :: dx2
real(wp), dimension(0:lx3+2) :: dx3
real(wp), dimension(1:lx1) :: dx1i
real(wp), dimension(1:lx2) :: dx2i
real(wp), dimension(1:lx3) :: dx3i

!! boundary condition arrays
real(wp), dimension(lx3) :: Vminx2,Vmaxx2
real(wp), dimension(lx2) :: Vminx3,Vmaxx3

!! mpi stuff
integer :: ierr,myid,lid

!! solution arrays
real(wp), allocatable, dimension(:,:) :: Phi

!! coefficient arrays
real(wp), allocatable, dimension(:,:) :: A, Ap, B, C, SigH

!! RHS array
real(wp), allocatable, dimension(:,:) :: srcterm
!real(wp), allocatable, dimension(:,:,:) :: srcterm2

!! MUMPS stuff
logical :: perflag=.false.     !shouldn't be used
integer :: it=1                !not used
real(wp) :: dt=1               !not used
integer :: gridflag=1
integer, dimension(4) :: flagsdirich=[1,1,1,1]        !denoting all Dirichlet conditions for test problem

!! command line input
character(4096) :: argv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! allocations
allocate(Phi(lx2,lx3))
allocate(A(lx2,lx3), Ap(lx2,lx3), SigH(lx2,lx3), B(lx2,lx3), C(lx2,lx3))
allocate(srcterm(lx2,lx3))

!! mpi stuff
call mpi_init()
call mpi_comm_rank(MPI_COMM_WORLD,myid)
call mpi_comm_size(MPI_COMM_WORLD,lid)

!! Set things up to give debug output
debug=.true.

!! Set up grid and compute differences needed for solution of PDE
call set_coordinates(x2dist, x3dist, &
                     x1,dx1,x1i,dx1i, &
                     x2,dx2,x2i,dx2i, &
                     x3,dx3,x3i,dx3i)

!! get coefficients and source terms
call set_parameters(x2,x3,A,Ap,SigH,B,C,srcterm,Vminx2,Vmaxx2,Vminx3,Vmaxx3,flagsdirich)

!! Make the call to GEMINI wrapper for PDE elliptic solver library, note the separate calls for root vs. workers
if (myid==0) then
  print*, 'Starting MUMPS solve...'
    Phi=elliptic2D_static(srcterm,A,Ap,SigH,B,C,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                                 dt,dx1,dx1i,dx2,dx2i,dx3,dx3i, &
                                 flagsdirich,perflag,it)
  print*, 'MUMPS solve is complete...'
else
  call elliptic_workers()
end if

!! Write some output for visualizations
if (myid==0) then
  call get_command_argument(1, argv, status=ierr)
  if(ierr /= 0) error stop 'please specify filename'

  print*, 'Numerical solution range:  ',minval(Phi),maxval(Phi)

  call hout%open(trim(argv), action="w")
  call hout%write("/lx1", lx1)
  call hout%write("/lx2", lx2)
  call hout%write("/lx3", lx3)
  call hout%write("/x1", x1(1:lx1))
  call hout%write("/x2", x2(1:lx2))
  call hout%write("/x3", x3(1:lx3))
  call hout%write("/Phi", Phi)
  call hout%write("/A",A)
  call hout%write("/Ap",Ap)
  call hout%write("/SigH",SigH)
  call hout%write("/B",B)
  call hout%write("/C",C)
  call hout%write("/srcterm",srcterm)
  call hout%close()
end if

call mpi_finalize()

contains
  !! populate coordinate arrays
  subroutine set_coordinates(x2dist,x3dist,x1,dx1,x1i,dx1i, &
                     x2,dx2,x2i,dx2i, &
                     x3,dx3,x3i,dx3i)
    real(wp) :: x2dist,x3dist
    real(wp), dimension(-1:) :: x1
    real(wp), dimension(-1:) :: x2
    real(wp), dimension(-1:) :: x3
    real(wp), dimension(1:) :: x1i
    real(wp), dimension(1:) :: x2i
    real(wp), dimension(1:) :: x3i
    real(wp), dimension(0:) :: dx1
    real(wp), dimension(0:) :: dx2
    real(wp), dimension(0:) :: dx3
    real(wp), dimension(1:) :: dx1i
    real(wp), dimension(1:) :: dx2i
    real(wp), dimension(1:) :: dx3i
    integer :: lx1,lx2,lx3,ix1,ix2,ix3

    lx1=size(x1)-4
    lx2=size(x2)-4
    lx3=size(x3)-4

    !! some basic bwd diffs and interface values, make sure domain is centered about zero in x2,3
    x1=[ (real(ix1-1,wp)/real(lx1-1,wp), ix1=-1,lx1+2) ]
    dx1=x1(0:lx1+2)-x1(-1:lx1+1)
    x1i(1:lx1+1)=0.5*(x1(0:lx1)+x1(1:lx1+1))
    dx1i=x1i(2:lx1+1)-x1i(1:lx1)
    
    x2=[ (real(ix2-1,wp)/real(lx2-1,wp), ix2=-1,lx2+2) ]*x2dist - x2dist/2.0
    dx2=x2(0:lx2+2)-x2(-1:lx2+1)
    x2i(1:lx2+1)=0.5*(x2(0:lx2)+x2(1:lx2+1))
    dx2i=x2i(2:lx2+1)-x2i(1:lx2)
    
    x3=[ (real(ix3-1,wp)/real(lx3-1,wp), ix3=-1,lx3+2) ]*x3dist - x3dist/2.0
    dx3=x3(0:lx3+2)-x3(-1:lx3+1)
    x3i(1:lx3+1)=0.5*(x3(0:lx3)+x3(1:lx3+1))
    dx3i=x3i(2:lx3+1)-x3i(1:lx3)
  end subroutine set_coordinates

  !! populate coefficient arrays
  subroutine set_parameters(x2,x3,A,Ap,SigH,B,C,srcterm,Vminx2,Vmaxx2,Vminx3,Vmaxx3,flagsdirich)
    real(wp), dimension(-1:) :: x2
    real(wp), dimension(-1:) :: x3
    real(wp), dimension(:,:) :: A, Ap, B, C, srcterm, SigH
    real(wp), dimension(:) :: Vminx2,Vmaxx2
    real(wp), dimension(:) :: Vminx3,Vmaxx3
    integer, dimension(4) :: flagsdirich
    integer :: lx2,lx3
    real(wp), parameter :: ell2=15e3
    real(wp), parameter :: ell3=100e3
    
    lx2=size(A,1); lx3=size(A,2);

    do ix3=1,lx3
      do ix2=1,lx2
        !! This is Pedersen conductance
        A(ix2,ix3)=0.1 + 10.0*(  exp(-(x2(ix2)-1.5*ell2)**2/2/ell2**2)*exp(-x3(ix3)**2/2/ell3**2) )
        Ap(ix2,ix3)=A(ix2,ix3)
        
        !! RHS
        srcterm(ix2,ix3)=-1e-6*( exp(-(x2(ix2)-1.5*ell2)**2/2/ell2**2)*exp(-x3(ix3)**2/2/ell3**2) - &
                exp(-(x2(ix2)+1.5*ell2)**2/2/ell2**2)*exp(-x3(ix3)**2/2/ell3**2) )
      end do
    end do
    SigH(1:lx2,1:lx3)=-3.0*A(1:lx2,1:lx3)
    call grad2D(x2,x3,SigH,C,B)

    Vminx2(1:lx3)=0.0; Vmaxx2(1:lx2)=0.0; Vminx3(1:lx2)=0.0; Vmaxx3(1:lx2)=0.0
    flagsdirich=[1,1,1,1]
  end subroutine set_parameters

  !! compute a simple spatial gradient (5 point stencil)
  subroutine grad2D(x2,x3,f,fx2,fx3)
    real(wp), dimension(-1:) :: x2
    real(wp), dimension(-1:) :: x3
    real(wp), dimension(:,:) :: f
    real(wp), dimension(:,:) :: fx2,fx3
    integer :: lx2,lx3,ix2,ix3

    lx2=size(x2)-4; lx3=size(x3)-4

    do ix3=1,lx3
      fx2(1,ix3) = (f(2,ix3)-f(1,ix3))/(x2(2)-x2(1))
      fx2(2:lx2-1,ix3) = (f(3:lx2,ix3)-f(1:lx2-2,ix3))/(x2(3:lx2)-x2(1:lx2-2))
      fx2(lx2,ix3) = (f(lx2,ix3)-f(lx2-1,ix3))/(x2(lx2)-x2(lx2-1))
    end do

    do ix2=1,lx2
      fx3(ix2,1) = (f(ix2,2)-f(ix2,1))/(x3(2)-x3(1))
      fx3(ix2,2:lx3-1) = (f(ix2,3:lx3)-f(ix2,1:lx3-2))/(x3(3:lx3)-x3(1:lx3-2))
      fx2(lx2,lx3) = (f(ix2,lx3)-f(ix2,lx3-1))/(x3(lx3)-x3(lx3-1))     
    end do
  end subroutine grad2D
end program test_potential2d_auroral

