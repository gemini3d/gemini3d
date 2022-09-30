submodule (io) plasma_input
!! plasma.f90 uses submodules in plasma_input_*.f90 and plasma_output_*.f90 for file I/O
use reader, only : get_simsize3
use sanity_check, only : check_finite_current, check_finite_plasma
use interpolation, only : interp3
use grid, only : get_grid3_coords_hdf5

implicit none (type, external)

interface ! plasma_input_*.f90

  module subroutine input_root_currents_hdf5(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
    character(*), intent(in) :: outdir
    integer, intent(in) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J1,J2,J3
    !! intent(out)
  end subroutine

  module subroutine input_root_mpi_hdf5(x1,x2all,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)
    real(wp), dimension(-1:), intent(in) :: x1, x2all, x3all
    character(*), intent(in) :: indatsize, indatfile
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,vs1,Ts
    !! intent(out)
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: Phi
    !! intent(out)
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: Phiall
    !! intent(out)
  end subroutine

  module subroutine getICs_hdf5(indatsize,indatfile,nsall,vs1all,Tsall,Phiall)
    character(*), intent(in) :: indatsize, indatfile
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: nsall,vs1all,Tsall
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: Phiall
  end subroutine

end interface

contains

subroutine input_root_currents(outdir,out_format, flagoutput,ymd,UTsec,J1,J2,J3)
  character(*), intent(in) :: outdir, out_format
  integer, intent(in) :: flagoutput
  integer, dimension(3), intent(in) :: ymd
  real(wp), intent(in) :: UTsec
  real(wp), dimension(:,:,:), intent(inout) :: J1,J2,J3
  !! intent(out)


  select case(out_format)
  case('h5')
    call input_root_currents_hdf5(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
  case default
    error stop 'input_root_current: unexpected Gemini input: ' // out_format
  end select
end subroutine input_root_currents


module procedure input_plasma
  ! subroutine input_plasma(x1,x2,x3all,indatsize,indatfile, ns,vs1,Ts)
  !! A BASIC WRAPPER FOR THE ROOT AND WORKER INPUT FUNCTIONS
  !! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
  !! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

  if (mpi_cfg%myid==0) then
    !! ROOT FINDS/CALCULATES INITIAL CONDITIONS AND SENDS TO WORKERS
    call input_root_mpi_hdf5(x1,x2,x3all,indatsize,indatfile,ns,vs1,Ts,Phi,Phiall)

    !> USER SUPPLIED FUNCTION TO TAKE A REFERENCE PROFILE AND CREATE INITIAL CONDITIONS FOR ENTIRE GRID.
    !> ASSUMING THAT THE INPUT DATA ARE EXACTLY THE CORRECT SIZE (AS IS THE CASE WITH FILE INPUT) THIS IS NOW SUPERFLUOUS
    print '(/,A,/,A)', 'Initial conditions (root):','------------------------'
    print '(A,2ES11.2)', 'Min/max input density:',     minval(ns(1:lx1,1:lx2,1:lx3,7)),  maxval(ns(1:lx1,1:lx2,1:lx3,7))
    print '(A,2ES11.2)', 'Min/max input velocity:',    minval(vs1(1:lx1,1:lx2,1:lx3,1:lsp)), maxval(vs1(1:lx1,1:lx2,1:lx3,1:lsp))
    print '(A,2ES11.2)', 'Min/max input temperature:', minval(Ts(1:lx1,1:lx2,1:lx3,1:lsp)),  maxval(Ts(1:lx1,1:lx2,1:lx3,1:lsp))
    print '(A,2ES11.2)', 'Min/max input electric potential:', minval(Phi(1:lx1,1:lx2,1:lx3)),  maxval(Phi(1:lx1,1:lx2,1:lx3))
    print '(A,2ES11.2)', 'Min/max input electric potential (full grid):', minval(Phiall(1:lx1,1:lx2all,1:lx3all)),  &
                         maxval(Phiall(1:lx1,1:lx2all,1:lx3all))

    call check_finite_plasma(out_dir, ns, vs1, Ts)
  else
    !! WORKERS RECEIVE THE IC DATA FROM ROOT
    call input_workers_mpi(ns,vs1,Ts,Phi)
  end if
end procedure input_plasma


!> Interpolate initial conditions onto "local" subgrid; we assume that the input data grid is specified
!    by the input file, whereas the target grid *could* be different, e.g. due to refinement or some other
!    custom arrangement.  The entire input file will be read by each worker calling this procedure.
!  Since this is only performing spatial interpolation it is easiest to just use the interpolation module
!    directly rather than create a type extension for inputdata (which inherently wants to also do time interpolation)
!    and then overriding the interp to space-only.
module procedure interp_file2subgrid
  real(wp), dimension(:,:,:,:), allocatable :: nsall,vs1all,Tsall
  real(wp), dimension(:,:,:), allocatable :: Phiall
  real(wp), dimension(:), allocatable :: parmflat
  integer :: lx1,lx2,lx3
  integer :: lx1in,lx2in,lx3in
  real(wp), dimension(:), allocatable :: x1in,x2in,x3in
  real(wp) :: glatctr=0._wp, glonctr=0._wp
  integer :: isp

  ! convenience
  lx1=size(x1)-4; lx2=size(x2)-4; lx3=size(x3)-4;

  ! read in the ICs size and allocate data
  call get_simsize3(out_dir // "/simsize.h5", lx1in,lx2in,lx3in)
  allocate(x1in(-1:lx1in+2),x2in(-1:lx2in+2),x3in(-1:lx3in+2))
  allocate(nsall(-1:lx1in+2,-1:lx2in+2,-1:lx3in+2,1:lsp), &
            vs1all(-1:lx1in+2,-1:lx2in+2,-1:lx3in+2,1:lsp), &
            Tsall(-1:lx1in+2,-1:lx2in+2,-1:lx3in+2,1:lsp), &
            Phiall(-1:lx1in+2,-1:lx2in+2,-1:lx3in+2))
  allocate(parmflat(lx1*lx2*lx3))

  ! get the input grid coordinates
  call get_grid3_coords_hdf5(out_dir,x1in,x2in,x3in,glonctr,glatctr)


  ! we must make sure that the target coordinates do not range outside the input file coordinates
  if(x1(1)<x1in(1) .or. x1(lx1)>x1in(lx1in)) then
    error stop 'interp_file2grid: x1 target coordinates beyond input grid coords'
  end if
  if(x2(1)<x2in(1) .or. x2(lx2)>x2in(lx2in)) then
    error stop 'interp_file2grid: x2 target coordinates beyond input grid coords'
  end if
  if(x3(1)<x3in(1) .or. x3(lx3)>x3in(lx3in)) then
    error stop 'interp_file2grid: x3 target coordinates beyond input grid coords'
  end if

  ! read in the input initial conditions, only hdf5 files are support for this functionality
  call getICs_hdf5(indatsize,indatfile,nsall,vs1all,Tsall,Phiall)

  ! interpolation input data to mesh sites; do not interpolate to ghost cells
  do isp=1,lsp
    parmflat=interp3(x1in(1:lx1in),x2in(1:lx2in),x3in(1:lx3in),nsall(1:lx1in,1:lx2in,1:lx3in,isp), &
                       x1(1:lx1),x2(1:lx2),x3(1:lx3))
    ns(1:lx1,1:lx2,1:lx3,isp)=reshape(parmflat,[lx1,lx2,lx3])
    parmflat=interp3(x1in(1:lx1in),x2in(1:lx2in),x3in(1:lx3in),vs1all(1:lx1in,1:lx2in,1:lx3in,isp), &
                       x1(1:lx1),x2(1:lx2),x3(1:lx3))
    vs1(1:lx1,1:lx2,1:lx3,isp)=reshape(parmflat,[lx1,lx2,lx3])
    parmflat=interp3(x1in(1:lx1in),x2in(1:lx2in),x3in(1:lx3in),Tsall(1:lx1in,1:lx2in,1:lx3in,isp), &
                       x1(1:lx1),x2(1:lx2),x3(1:lx3))
    Ts(1:lx1,1:lx2,1:lx3,isp)=reshape(parmflat,[lx1,lx2,lx3])
  end do
  parmflat=interp3(x1in(1:lx1in),x2in(1:lx2in),x3in(1:lx3in),Phiall(1:lx1in,1:lx2in,1:lx3in), &
                     x1(1:lx1),x2(1:lx2),x3(1:lx3))
  Phi(1:lx1,1:lx2,1:lx3)=reshape(parmflat,[lx1,lx2,lx3])

  deallocate(x1in,x2in,x3in,nsall,vs1all,Tsall,Phiall,parmflat)
end procedure interp_file2subgrid


module procedure input_plasma_currents
  ! module subroutine input_plasma_currents(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
  !! READS, AS INPUT, A FILE GENERATED BY THE GEMINI.F90 PROGRAM.
  !! THIS SUBROUTINE IS A WRAPPER FOR SEPARATE ROOT/WORKER CALLS

  if (mpi_cfg%myid==0) then
    !> ROOT FINDS/CALCULATES INITIAL CONDITIONS AND SENDS TO WORKERS
    print *, 'Assembling current density data on root...  '
    call input_root_currents(outdir, out_format,flagoutput,ymd,UTsec,J1,J2,J3)

    call check_finite_current(outdir, J1, J2, J3)
  else
    !> WORKERS RECEIVE THE IC DATA FROM ROOT
    call input_workers_currents(J1,J2,J3)
  end if
end procedure input_plasma_currents


subroutine input_workers_currents(J1,J2,J3)
  !! WORKER INPUT FUNCTIONS FOR GETTING CURRENT DENSITIES

  real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J1,J2,J3
  !! intent(out)

  !> ALL WE HAVE TO DO IS WAIT TO RECEIVE OUR PIECE OF DATA FROM ROOT
  call bcast_recv3D_ghost(J1,tag%J1)
  call bcast_recv3D_ghost(J2,tag%J2)
  call bcast_recv3D_ghost(J3,tag%J3)
end subroutine input_workers_currents


subroutine input_workers_mpi(ns,vs1,Ts,Phi)
  !------------------------------------------------------------
  !-------RECEIVE INITIAL CONDITIONS FROM ROOT PROCESS
  !------------------------------------------------------------

  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,vs1,Ts
  !! intent(out)

  real(wp), dimension(-1:,-1:,-1:) :: Phi

  ns=0
  vs1=0
  Ts=0

  call bcast_recv(ns,tag%ns)
  call bcast_recv(vs1,tag%vs1)
  call bcast_recv(Ts,tag%Ts)
  call bcast_recv3D_ghost(Phi,tag%Phi)

  ! print*, mpi_cfg%myid
  ! print *, 'Min/max input density:  ',     minval(ns(:,:,:,7)),  maxval(ns(:,:,:,7))
  ! print *, 'Min/max input velocity:  ',    minval(vs1(:,:,:,:)), maxval(vs1(:,:,:,:))
  ! print *, 'Min/max input temperature:  ', minval(Ts(:,:,:,:)),  maxval(Ts(:,:,:,:))
end subroutine input_workers_mpi


end submodule plasma_input
