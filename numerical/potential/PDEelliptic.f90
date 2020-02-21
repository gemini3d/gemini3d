module PDEelliptic

!! Various tools for solving elliptic partial differential equations - uses MUMPS, scalapack, lapack, openmpi, and blas

use, intrinsic:: iso_fortran_env, only: stderr=>error_unit, stdout=>output_unit
use mumps_interface, only : mumps_struc, mumps_exec
use mpi, only: mpi_comm_world
use phys_consts, only: wp, debug

implicit none
private
public :: elliptic3D_cart,elliptic2D_cart,elliptic2D_polarization,elliptic2D_polarization_periodic,&
  elliptic_workers, check_mumps_status, quiet_mumps

interface ! elliptic2d.f90

module function elliptic2D_polarization(srcterm,SigP2,SigP3,SigH,gradSigH2,gradSigH3,Cm,v2,v3,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dt,dx1, &
  dx1i,dx2all,dx2iall,dx3all,dx3iall,Phi0,perflag,it)
real(wp), dimension(:,:), intent(in) :: srcterm,SigP2,SigP3,SigH,Cm,gradSigH2,gradSigH3,v2,v3
!! ZZZ - THESE WILL NEED TO BE MODIFIED CONDUCTIVITIES, AND WE'LL NEED THREE OF THEM
real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
real(wp), dimension(:), intent(in) :: Vminx3,Vmaxx3
real(wp), intent(in) :: dt
real(wp), dimension(0:), intent(in) :: dx1         !backward diffs start at index zero due to ghost cells
real(wp), dimension(:), intent(in) :: dx1i         !centered diffs do not include any ghost cells
real(wp), dimension(0:), intent(in) :: dx2all
real(wp), dimension(:), intent(in) :: dx2iall
real(wp), dimension(0:), intent(in) :: dx3all
real(wp), dimension(:), intent(in) :: dx3iall
real(wp), dimension(:,:), intent(in) :: Phi0
logical, intent(in) :: perflag
integer, intent(in) :: it
real(wp), dimension(size(SigP2,1),size(SigP2,2)) :: elliptic2D_polarization
end function elliptic2D_polarization

module function elliptic2D_polarization_periodic(srcterm,SigP,SigH,gradSigH2,gradSigH3,Cm,v2,v3,Vminx2,Vmaxx2, &
    Vminx3,Vmaxx3,dt,dx1,dx1i,dx2all,dx2iall,dx3all,dx3iall,Phi0,perflag,it)
real(wp), dimension(:,:), intent(in) :: srcterm,SigP,SigH,gradSigH2,gradSigH3,Cm,v2,v3
real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
real(wp), dimension(:), intent(in) :: Vminx3,Vmaxx3
real(wp), intent(in) :: dt
real(wp), dimension(0:), intent(in) :: dx1         !backward diffs start at index zero due to ghost cells
real(wp), dimension(:), intent(in) :: dx1i         !centered diffs do not include any ghost cells
real(wp), dimension(0:), intent(in) :: dx2all
real(wp), dimension(:), intent(in) :: dx2iall
real(wp), dimension(0:), intent(in) :: dx3all
real(wp), dimension(:), intent(in) :: dx3iall
real(wp), dimension(:,:), intent(in) :: Phi0
logical, intent(in) :: perflag
integer, intent(in) :: it
real(wp), dimension(size(SigP,1),size(SigP,2)) :: elliptic2D_polarization_periodic
end function elliptic2D_polarization_periodic

module function elliptic2D_cart(srcterm,sig0,sigP,Vminx1,Vmaxx1,Vminx3,Vmaxx3,&
  dx1,dx1i,dx3all,dx3iall,flagdirich,perflag,gridflag,it)
real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP   !arrays passed in will still have full rank 3
real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
real(wp), dimension(0:), intent(in) :: dx1         !backward diffs start at index zero due to ghost cells
real(wp), dimension(:), intent(in) :: dx1i         !centered diffs do not include any ghost cells
real(wp), dimension(0:), intent(in) :: dx3all
real(wp), dimension(:), intent(in) :: dx3iall
integer, intent(in) :: flagdirich
logical, intent(in) :: perflag
integer, intent(in) :: gridflag
integer, intent(in) :: it
real(wp), dimension(size(sig0,1),1,size(sig0,3)) :: elliptic2D_cart
end function elliptic2D_cart

end interface


integer, dimension(:), pointer, protected, save :: mumps_perm   !cached permutation, unclear whether save is necessary...




contains


function elliptic3D_cart(srcterm,Ac,Bc,Cc,Dc,Ec,Fc,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                  dx1,dx1i,dx2all,dx2iall,dx3all,dx3iall,flagdirich,perflag,it)

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
!-------ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
!-------LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
!-------GRID POINTS.  This version requires that one compute the
!-------Coefficients of the problem beforehand and solves the equation:
!-------
!-------A d^2V/dx1^2 + B d^2V/dx2^2 + C d^2V/dx3^2 + D dV/dx1 + E dV/dx2 + F dV/dx3 = srcterm
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

type (MUMPS_STRUC) :: mumps_par

real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D_cart

!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
!if (myid==0) then
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

  if (debug) print *, 'MUMPS will attempt a solve of size:  ',lx1,lx2,lx3
  if (debug) print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


  !! DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER GUIDE).
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
!end if
if (debug) print *, 'Number of entries used:  ',ient-1


! INIT MUMPS

mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call MUMPS_exec(mumps_par)

call quiet_mumps(mumps_par)

!LOAD OUR PROBLEM
if (debug) print*, 'Loading mumps problem...'
!if ( myid==0 ) then
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

  if (debug) print*, 'Dealing with permutation',perflag,it
  if (perflag .and. it/=1) then       !used cached permutation, but only gets tested on first time step
    allocate(mumps_par%PERM_IN(mumps_par%N))
    mumps_par%PERM_IN=mumps_perm
    mumps_par%ICNTL(7)=1
  end if


  if (debug) print*, 'Setting memory relaxation...'
  !3D solves very often need better memory relaxation
  mumps_par%ICNTL(14)=500
!end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
if (debug) print*, 'Executing solve...'
mumps_par%JOB = 6

call MUMPS_exec(mumps_par)

call check_mumps_status(mumps_par, 'elliptic3D_cart')

!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
!if ( myid==0 ) then
  if (debug) print *, 'Now organizing results...'

  if (perflag .and. it==1) then
    allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
    mumps_perm=mumps_par%SYM_PERM
  end if

  elliptic3D_cart=reshape(mumps_par%RHS,[lx1,lx2,lx3])

  if (debug) print *, 'Now attempting deallocations...'

  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
!end if

mumps_par%JOB = -2

call MUMPS_exec(mumps_par)

end function elliptic3D_cart




subroutine quiet_mumps(obj)
!! this must be called AFTER the first mumps call that had job=-1
!! Needs Mumps >= 5.2 to actually take effect, does nothing on older MUMPS
!! it stops the 100's of megabytes of logging console text, probably speeding up as well

type(MUMPS_STRUC), intent(inout) :: obj

obj%icntl(1) = stderr  ! error messages
obj%icntl(2) = stdout !  diagnosic, statistics, and warning messages
obj%icntl(3) = stdout! ! global info, for the host (myid==0)
obj%icntl(4) = 1           ! default is 2, this reduces verbosity

end subroutine quiet_mumps


subroutine elliptic_workers()
!! ALLOWS WORKERS TO ENTER MUMPS SOLVES

type(MUMPS_STRUC) :: mumps_par

!FIRE UP MUMPS
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call MUMPS_exec(mumps_par)

call quiet_mumps(mumps_par)


!ROOT WILL LOAD OUR PROBLEM


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6

call MUMPS_exec(mumps_par)

call check_mumps_status(mumps_par, 'elliptic_workers')

!DEALLOCATE STRUCTURES USED BY WORKERS DURING SOLVE
mumps_par%JOB = -2

call MUMPS_exec(mumps_par)

end subroutine elliptic_workers


subroutine check_mumps_status(p, name)
!! check if Mumps error occurred

type(MUMPS_STRUC), intent(in) :: p
character(*), intent(in) :: name

if (p%info(1) < 0 .or. p%infog(1) < 0) then
  write(stderr, *) 'Gemini:PDEelliptic:' // name // '  MUMPS ERROR: details:'
  if (p%info(1) == -1) write(stderr,'(a,i4)') 'the error was reported by processor #',p%info(2)
  write(stderr, *) 'Mumps Error: info(1,2,8):', p%info(1:2), p%info(8)
  write(stderr, *) 'Mumps Error: infoG(1,2)', p%infoG(1:2)
  write(stderr, *) 'for error number meaning, see "8 Error Diagnostics" of MUMPS manual'
  error stop
endif

end subroutine check_mumps_status

end module PDEelliptic
