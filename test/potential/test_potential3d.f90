program test_potential3D

use, intrinsic :: iso_fortran_env, only: real64

use mpi_f08, only : mpi_init, mpi_comm_rank, MPI_COMM_WORLD,mpi_finalize
use mumps_interface, only : mumps_struc, mumps_exec
use h5fortran, only : hdf5_file

implicit none (type, external)

type(mumps_struc) :: mumps_par

type(hdf5_file) :: hout

integer :: ierr, myid

! integer, parameter :: npts1=34,npts2=144,npts3=96
integer, parameter :: npts1=14, npts2=24, npts3=16
!! made smaller for unit tests

integer, parameter :: lk=npts1*npts2*npts3
integer :: lent
integer :: ix1,ix2,ix3, lx1,lx2,lx3
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(real64), dimension(:), allocatable :: M
real(real64), dimension(:), allocatable :: b
real(real64) :: dx1
real(real64), dimension(npts2,npts3) :: Vminx1,Vmaxx1
real(real64), dimension(npts1,npts3) :: Vminx2,Vmaxx2
real(real64), dimension(npts1,npts2) :: Vminx3,Vmaxx3
real(real64), dimension(:,:), allocatable ::  Mfull
real(real64) :: tstart,tfin

character(4096) :: argv

call MPI_INIT()

call mpi_comm_rank(MPI_COMM_WORLD,myid)

!------------------------------------------------------------
!-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
!-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
!-------GUIDE).
!------------------------------------------------------------
lent=7*(npts1-2)*(npts2-2)*(npts3-2)                                           !interior entries
lent=lent+2*(npts1-2)*(npts2-2)+2*(npts2-2)*(npts3-2)+2*(npts1-2)*(npts3-2)    !6 faces of cube
lent=lent+4*(npts1-2)+4*(npts2-2)+4*(npts3-2)                                  !12 edges
lent=lent+8                                                                    !8 corners
allocate(ir(lent),ic(lent),M(lent),b(lk))
lx1=npts1
lx2=npts2
lx3=npts3

dx1 = 1. / npts1           !scale dx so the domain of problem is [0,1]

Vminx1(:,:) = 0
Vmaxx1(:,:) = 0
Vminx2(:,:) = 0
Vmaxx2(:,:) = 0
Vminx3(:,:) = 0
Vmaxx3(:,:) = 10

M(:) = 0
b(:) = 0
ient=1


!LOAD UP MATRIX ELEMENTS
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      iPhi=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix2,ix3) as a column vector.  Also row # of big matrix

      if (ix1==1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient) = 1
        b(iPhi)=Vminx1(ix2,ix3)
        ient=ient+1
      elseif (ix1==lx1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient) = 1
        b(iPhi)=Vmaxx1(ix2,ix3)
        ient=ient+1
      elseif (ix2==1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient) = 1
        b(iPhi)=Vminx2(ix1,ix3)
        ient=ient+1
      elseif (ix2==lx2) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient) = 1
        b(iPhi)=Vmaxx2(ix1,ix3)
        ient=ient+1
      elseif (ix3==1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient) = 1
        b(iPhi)=Vminx3(ix1,ix2)
        ient=ient+1
      elseif (ix3==lx3) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient) = 1
        b(iPhi)=Vmaxx3(ix1,ix2)
        ient=ient+1
      else                       !INTERIOR
        !ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1*lx2
        M(ient) = 1
        ient=ient+1

        !ix1,ix2-1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1
        M(ient) = 1
        ient=ient+1

        !ix1-1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-1
        M(ient) = 1
        ient=ient+1

        !ix1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient) = -6
        ient=ient+1

        !ix1+1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient) = 1
        ient=ient+1

        !ix1,ix2+1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1
        M(ient) = 1
        ient=ient+1

        !ix1,ix2,ix3+1
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1*lx2
        M(ient) = 1
        ient=ient+1
      end if
    end do
  end do
end do


!CORRECT FOR DX /= 1
b=b*dx1**2


!OUTPUT FULL MATRIX FOR DEBUGGING IF ITS NOT TOO BIG (ZZZ --> CAN BE COMMENTED OUT)
if (myid==0) then
  call get_command_argument(1, argv, status=ierr)
  if(ierr /= 0) error stop 'please specify filename'

  call hout%open(trim(argv), action="w")
  call hout%write("/lx1", lx1)
  call hout%write("/lx2", lx2)
  call hout%write("/lx3", lx3)
  call hout%write("/b", b)
end if


if (lk<150) then
  allocate(Mfull(lk,lk))
  Mfull(:,:) = 0
  do ient=1,size(ir)
    Mfull(ir(ient),ic(ient))=M(ient)
  end do

  if (myid==0) call hout%write("/Mfull", Mfull)

  deallocate(Mfull)
end if


!------------------------------------------------------------
!-------DO SOME STUFF TO CALL MUMPS
!------------------------------------------------------------

! Define a communicator for the package.
mumps_par%COMM = MPI_COMM_WORLD%mpi_val


!Initialize an instance of the package
!for L U factorization (sym = 0, with working host)
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
call mumps_exec(mumps_par)


!Define problem on the host (processor 0)
if ( mumps_par%MYID == 0 ) then
  mumps_par%N=lk
  mumps_par%NZ=lent
  allocate( mumps_par%IRN ( mumps_par%NZ ) )
  allocate( mumps_par%JCN ( mumps_par%NZ ) )
  allocate( mumps_par%A( mumps_par%NZ ) )
  allocate( mumps_par%RHS ( mumps_par%N  ) )
  mumps_par%IRN=ir
  mumps_par%JCN=ic
  mumps_par%A=M
  mumps_par%RHS=b

!  mumps_par%ICNTL(7)=6    !force a particular reordering - see mumps docs
!  mumps_par%ICNTL(28)=2
!  mumps_par%ICNTL(29)=2
end if


!Call package for solution
mumps_par%JOB = 6
call cpu_time(tstart)
call mumps_exec(mumps_par)
call cpu_time(tfin)
write(*,*) 'Solve took ',tfin-tstart,' seconds...'


!Solution has been assembled on the host
if ( mumps_par%MYID == 0 ) then
  call hout%write("/x1", mumps_par%RHS/dx1**2)
end if

call hout%close()

!Deallocate user data
if ( mumps_par%MYID == 0 ) then
  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
end if
deallocate(ir,ic,M,b)


!Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
call mumps_exec(mumps_par)

call MPI_FINALIZE()

end program
