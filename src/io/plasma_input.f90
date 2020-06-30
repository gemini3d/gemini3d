submodule (io) plasma_input
!! plasma.f90 uses submodules in plasma_input_*.f90 and plasma_output_*.f90 for raw, hdf5 or netcdf4 I/O
use reader, only : get_simsize3
use pathlib, only : get_suffix
use sanity_check, only : check_finite_current, check_finite_plasma

implicit none (type, external)

interface ! plasma_input_*.f90

module subroutine input_root_currents_raw(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
character(*), intent(in) :: outdir
integer, intent(in) :: flagoutput
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(out) :: J1,J2,J3
end subroutine input_root_currents_raw

module subroutine input_root_mpi_raw(x1,x2all,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)
real(wp), dimension(-1:), intent(in) :: x1, x2all, x3all
character(*), intent(in) :: indatsize, indatfile
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts
real(wp), dimension(:,:,:), intent(out) :: Phi
real(wp), dimension(:,:,:), intent(out) :: Phiall
end subroutine input_root_mpi_raw

module subroutine input_root_currents_hdf5(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
character(*), intent(in) :: outdir
integer, intent(in) :: flagoutput
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(out) :: J1,J2,J3
end subroutine input_root_currents_hdf5

module subroutine input_root_mpi_hdf5(x1,x2all,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)
real(wp), dimension(-1:), intent(in) :: x1, x2all, x3all
character(*), intent(in) :: indatsize, indatfile
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts
real(wp), dimension(:,:,:), intent(out) :: Phi
real(wp), dimension(:,:,:), intent(out) :: Phiall
end subroutine input_root_mpi_hdf5

module subroutine input_root_currents_nc4(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
character(*), intent(in) :: outdir
integer, intent(in) :: flagoutput
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(out) :: J1,J2,J3
end subroutine input_root_currents_nc4

module subroutine input_root_mpi_nc4(x1,x2all,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)
real(wp), dimension(-1:), intent(in) :: x1, x2all, x3all
character(*), intent(in) :: indatsize, indatfile
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts
real(wp), dimension(:,:,:), intent(out) :: Phi
real(wp), dimension(:,:,:), intent(out) :: Phiall
end subroutine input_root_mpi_nc4

end interface

contains

subroutine input_root_currents(outdir,out_format, flagoutput,ymd,UTsec,J1,J2,J3)
character(*), intent(in) :: outdir, out_format
integer, intent(in) :: flagoutput
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(out) :: J1,J2,J3


select case(out_format)
case('dat')
  call input_root_currents_raw(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
case('h5')
  call input_root_currents_hdf5(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
case ('nc')
  call input_root_currents_nc4(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
case default
  error stop 'input_root_current: unexpected Gemini input'
end select

end subroutine input_root_currents


module procedure input_plasma
! subroutine input_plasma(x1,x2,x3all,indatsize,indatfile, ns,vs1,Ts)
!! A BASIC WRAPPER FOR THE ROOT AND WORKER INPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

if (myid==0) then
  !! ROOT FINDS/CALCULATES INITIAL CONDITIONS AND SENDS TO WORKERS
  select case (get_suffix(indatsize))
  case ('.h5')
    call input_root_mpi_hdf5(x1,x2,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)
  case ('.nc')   !neither netcdf now raw input support restarting right now
    call input_root_mpi_nc4(x1,x2,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)
  case ('.dat')
    call input_root_mpi_raw(x1,x2,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)
  case default
    write(stderr,*) 'input_plasma: unknown grid format: ' // get_suffix(indatsize)
    error stop 6
  end select

  !> USER SUPPLIED FUNCTION TO TAKE A REFERENCE PROFILE AND CREATE INITIAL CONDITIONS FOR ENTIRE GRID.
  !> ASSUMING THAT THE INPUT DATA ARE EXACTLY THE CORRECT SIZE (AS IS THE CASE WITH FILE INPUT) THIS IS NOW SUPERFLUOUS
  print '(/,A,/,A)', 'Initial conditions (root):','------------------------'
  print '(A,2ES11.2)', 'Min/max input density:',     minval(ns(:,:,:,7)),  maxval(ns(:,:,:,7))
  print '(A,2ES11.2)', 'Min/max input velocity:',    minval(vs1(:,:,:,:)), maxval(vs1(:,:,:,:))
  print '(A,2ES11.2)', 'Min/max input temperature:', minval(Ts(:,:,:,:)),  maxval(Ts(:,:,:,:))
  print '(A,2ES11.2)', 'Min/max input electric potential:', minval(Phi(:,:,:)),  maxval(Phi(:,:,:))
  print '(A,2ES11.2)', 'Min/max input electric potential (full grid):', minval(Phiall(:,:,:)),  maxval(Phiall(:,:,:))

  call check_finite_plasma(ns, vs1, Ts)
else
  !! WORKERS RECEIVE THE IC DATA FROM ROOT
  call input_workers_mpi(ns,vs1,Ts,Phi)
end if

end procedure input_plasma


module procedure input_plasma_currents
! module subroutine input_plasma_currents(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
!! READS, AS INPUT, A FILE GENERATED BY THE GEMINI.F90 PROGRAM.
!! THIS SUBROUTINE IS A WRAPPER FOR SEPARATE ROOT/WORKER CALLS

if (myid==0) then
  !> ROOT FINDS/CALCULATES INITIAL CONDITIONS AND SENDS TO WORKERS
  print *, 'Assembling current density data on root...  '
  call input_root_currents(outdir,out_format,flagoutput,ymd,UTsec,J1,J2,J3)

  call check_finite_current(J1, J2, J3)
else
  !> WORKERS RECEIVE THE IC DATA FROM ROOT
  call input_workers_currents(J1,J2,J3)
end if

end procedure input_plasma_currents


subroutine input_workers_currents(J1,J2,J3)
!! WORKER INPUT FUNCTIONS FOR GETTING CURRENT DENSITIES

real(wp), dimension(:,:,:), intent(out) :: J1,J2,J3


!> ALL WE HAVE TO DO IS WAIT TO RECEIVE OUR PIECE OF DATA FROM ROOT
call bcast_recv(J1,tag%J1)
call bcast_recv(J2,tag%J2)
call bcast_recv(J3,tag%J3)

end subroutine input_workers_currents

subroutine input_workers_mpi(ns,vs1,Ts,Phi)

!------------------------------------------------------------
!-------RECEIVE INITIAL CONDITIONS FROM ROOT PROCESS
!------------------------------------------------------------

real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts
real(wp), dimension(:,:,:) :: Phi

ns=0
vs1=0
Ts=0

call bcast_recv(ns,tag%ns)
call bcast_recv(vs1,tag%vs1)
call bcast_recv(Ts,tag%Ts)
call bcast_recv(Phi,tag%Phi)

if (.false.) then
  print*, myid
  print *, 'Min/max input density:  ',     minval(ns(:,:,:,7)),  maxval(ns(:,:,:,7))
  print *, 'Min/max input velocity:  ',    minval(vs1(:,:,:,:)), maxval(vs1(:,:,:,:))
  print *, 'Min/max input temperature:  ', minval(Ts(:,:,:,:)),  maxval(Ts(:,:,:,:))
endif

end subroutine input_workers_mpi


end submodule plasma_input
