module PDEelliptic

!! Various tools for solving elliptic partial differential equations - uses MUMPS, scalapack, lapack, openmpi, and blas

use, intrinsic:: iso_fortran_env, only: stderr=>error_unit, stdout=>output_unit
use mumps_interface, only : mumps_struc, mumps_exec
use phys_consts, only: wp, debug

use mpi_f08, only: MPI_COMM_WORLD

implicit none (type, external)
private
public :: elliptic3D_cart,elliptic2D_cart,elliptic2D_polarization,elliptic2D_polarization_periodic,&
  elliptic_workers, check_mumps_status, quiet_mumps, elliptic3D_cart_periodic, &
  mumps_perm

interface ! elliptic2d.f90
  module function elliptic2D_polarization(srcterm,SigP2,SigP3,SigH,gradSigH2,gradSigH3,Cm,v2,v3,Vminx2,Vmaxx2, &
      Vminx3,Vmaxx3,dt,dx1, &
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
      dx1,dx1i,dx3all,dx3iall,flagsdirich,perflag,gridflag,it)
    real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP   !arrays passed in will still have full rank 3
    real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
    real(wp), dimension(0:), intent(in) :: dx1         !backward diffs start at index zero due to ghost cells
    real(wp), dimension(:), intent(in) :: dx1i         !centered diffs do not include any ghost cells
    real(wp), dimension(0:), intent(in) :: dx3all
    real(wp), dimension(:), intent(in) :: dx3iall
    integer, dimension(4), intent(in) :: flagsdirich      ! note this is used to set all of the boundary conditions for this specific routine:  [x1min,x1max,x3min,x3max]
    logical, intent(in) :: perflag
    integer, intent(in) :: gridflag
    integer, intent(in) :: it
    real(wp), dimension(size(sig0,1),size(sig0,2),size(sig0,3)) :: elliptic2D_cart
  end function elliptic2D_cart
end interface

interface ! elliptic3d
  module function elliptic3D_cart(srcterm,Ac,Bc,Cc,Dc,Ec,Fc,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
    dx1,dx1i,dx2all,dx2iall,dx3all,dx3iall,flagdirich,perflag,it)
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
    real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D_cart
  end function elliptic3D_cart

  module function elliptic3D_cart_periodic(srcterm,Ac,Bc,Cc,Dc,Ec,Fc,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
    dx1,dx1i,dx2all,dx2iall,dx3all,dx3iall,flagdirich,perflag,it)
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
    real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D_cart_periodic
  end function elliptic3D_cart_periodic
end interface

integer, dimension(:), pointer, protected :: mumps_perm

contains


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
  mumps_par%COMM = MPI_COMM_WORLD%mpi_val
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
