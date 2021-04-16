module meshobj

!> This is a base type for defining most general characteristics of a curvlinear mesh; extension define specific coordinate systems.
!    However, the idea here is to insulate the numerical parts of the program from that so that they can deal exclusively with the
!    generic metric factors, etc.

use phys_consts, only : wp
use h5fortran, only : hdf5_file
implicit none (type, external)
public

!> curvmesh is a "top-level" derived type containing functionality and data that is not specific to individual coordinate systems
!   (which are extended types).  Note that all arrays are pointers because they need to be targets and allocatable AND the fortran
!   standard does not support having allocatable, target attributed inside a derived type.  Because of this is it not straightforward
!   to check the allocation status of these arrays (i.e. fortran also does not allow one to check the allocation status of a pointer). 
!   Thus the quantities which are not, for sure, allocated need to have an allocation status variable so we can check...  Because
!   the pointers are always allocated in groups we do not need separate status vars for each array thankfully...
type :: curvmesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Generic properties !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !> we need to know whether or not different groups of pointers are allocated and the intrinsic "allocatable" only work on allocatables (not pointers)...
  logical :: xi_alloc_status=.false.
  logical :: dxi_alloc_status=.false. 
  logical :: difflen_alloc_status=.false.
  logical :: null_alloc_status=.false.
  logical :: geog_set_status=.false.      ! a distinction here is that geographic coords. get allocated with other arrays, but get set separately

  !> sizes.  Specified and set by base class methods
  integer :: lx1,lx2,lx3,lx2all,lx3all
  !! for program units that may not be able to access module globals
  
  !!> curvilinear vars. and diffs.
  real(wp), dimension(:), pointer :: x1     ! provided in input file
  real(wp), dimension(:), pointer :: x1i    ! recomputed by base class once x1 set
  real(wp), dimension(:), pointer :: dx1        ! recomputed by base class from x1
  !! because sub arrays need to be assigned to aliases in calculus module program units
  real(wp), dimension(:), pointer :: dx1i   ! recomputed by base class from interface locations (themselves recomputed)
  
  real(wp), dimension(:), pointer :: x2
  real(wp), dimension(:), pointer :: x2i
  real(wp), dimension(:), pointer :: dx2        ! this (and similar dx arrays) are pointers because they become associated with other pointers (meaning they have to either have the "target" keyword or themselves be pointers).  These should also be contiguous but I believe that is guaranteed as long as they are assigned through allocate statements (see fortran standard)  derived type arrays cannot be targets so we are forced to declare them as pointers and trust that they are allocated contiguous
  real(wp), dimension(:), pointer :: dx2i
  
  real(wp), dimension(:), pointer :: x3
  real(wp), dimension(:), pointer :: x3i
  real(wp), dimension(:), pointer :: dx3
  real(wp), dimension(:), pointer :: dx3i
  
  real(wp), dimension(:), pointer :: x2all
  real(wp), dimension(:), pointer  :: x2iall
  real(wp), dimension(:), pointer :: dx2all
  real(wp), dimension(:), pointer  :: dx2iall
  
  real(wp), dimension(:), pointer :: x3all
  real(wp), dimension(:), pointer  :: x3iall
  real(wp), dimension(:), pointer :: dx3all
  real(wp), dimension(:), pointer  :: dx3iall
  
  !> differential length elements.  Compute by a method in generic class given metric factors
  real(wp), dimension(:,:,:), pointer :: dl1i,dl2i,dl3i
  
  !> flag for indicating whether or not the grid is periodic
  logical :: flagper
  
  !> flag for indicated type of grid (0 - closed dipole; 1 - open dipole inverted; 2 - non-inverted).  Computed by method in generic
  !class once coordinate specific quantitaties are defined
  integer :: gridflag
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Coordinate system specific properties !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !> Metric factors.  These are pointers to be assigned/allocated/filled by type extensions for specific coordinate system
  logical :: coord_alloc_status=.false.    ! single status variable for all coord-specific arrays

  real(wp), dimension(:,:,:), pointer :: h1,h2,h3     ! need to be computed by subclass for specific coordinate system
  !! these are the cell-centered metric coefficients
  real(wp), dimension(:,:,:), pointer :: h1x1i,h2x1i,h3x1i
  !! metric factors at x1 cell interfaces; dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), pointer :: h1x2i,h2x2i,h3x2i
  !! metric factors at x2 interfaces; dim. 2 has lx2+1
  real(wp), dimension(:,:,:), pointer :: h1x3i,h2x3i,h3x3i
  !! metric factors at x3 interfaces; dim. 3 has lx3+1
  
  !> root only; full grid metric factors.  These must have type-extension bound procedures to provide allocation and filling
  real(wp), dimension(:,:,:), pointer :: h1all,h2all,h3all
  real(wp), dimension(:,:,:), pointer :: h1x1iall,h2x1iall,h3x1iall
  !! dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), pointer :: h1x2iall,h2x2iall,h3x2iall
  !! dim. 2 has lx2all+1
  real(wp), dimension(:,:,:), pointer :: h1x3iall,h2x3iall,h3x3iall
  !! dim. 3 has lx3all+1
  
  !> Problem-depedent geometric terms that can be precomputed, e.g. for advection and elliptic equations
  !! to be added in later if we deem it worthwhile to include these to save computation time and/or memory

  !> unit vectors - Pointers.  Subclass methods must allocate and assign values to these.  
  real(wp), dimension(:,:,:,:), pointer :: e1,e2,e3
  !! unit vectors in curvilinear space (cartesian components)
  real(wp), dimension(:,:,:,:), pointer :: er,etheta,ephi
  !! spherical unit vectors (cartesian components)
  
  !> geomagnetic grid data, ECEF spherical referenced to dipole axis
  real(wp), dimension(:,:,:), pointer :: r,theta,phi
  !! may be used in the interpolation of neutral perturbations
  
  !> root only; full grid geomagnetic positions.  Pointers.  Need type-extension bound procedures for each
  !coordinate system must allocate and compute these. 
  real(wp), dimension(:,:,:), pointer :: rall,thetaall,phiall
  
  !> geographic data; pointers
  real(wp), dimension(:,:,:), pointer :: glat,glon,alt
  
  !> magnetic field magnitude and inclination.  Pointers
  real(wp), dimension(:,:,:), pointer :: Bmag
  real(wp), dimension(:,:), pointer :: I

  !> Gravitational field
  real(wp), dimension(:,:,:), pointer :: g1,g2,g3
  
  !> EIA-required fullgrid data.  Pointers
  real(wp), dimension(:,:,:), pointer :: altall,glonall
  real(wp), dimension(:,:,:), pointer :: Bmagall
  
  !> points that should not be considered part of the numerical domain. Pointers
  logical, dimension(:,:,:), pointer :: nullpts
  !! this could be a logical but I'm going to treat it as real*8
  integer :: lnull
  !! length of null point index array
  integer, dimension(:,:), pointer :: inull

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! type-bound procedures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  contains
    procedure :: set_coords
    procedure :: calc_coord_diffs
    procedure :: calc_difflengths
    procedure :: calc_inull
    procedure :: calc_gridflag
    procedure :: init_storage
    procedure :: writesize
    procedure :: writegrid
    procedure :: writegridall
    final :: destructor
end type curvmesh

contains
  !> assign coordinates to internal variables given some set of input arrays.
  !   Assume that the data passed in include ghost cells
  subroutine set_coords(self,x1,x2,x3,x2all,x3all)
    class(curvmesh) :: self
    real(wp), dimension(:), intent(in) :: x1,x2,x3
    real(wp), dimension(:), intent(in) :: x2all,x3all
    integer :: lx1,lx2,lx3,lx2all,lx3all

    lx1=size(x1,1)-4
    lx2=size(x2,1)-4
    lx3=size(x3,1)-4
    lx2all=size(x2all,1)-4
    lx3all=size(x3all,1)-4
    self%lx1=lx1; self%lx2=lx2; self%lx3=lx3
    self%lx2all=lx2all; self%lx3all=lx3all

    allocate(self%x1(-1:lx1+2),self%x2(-1:lx2+2),self%x3(-1:lx3+2))
    self%x1=x1; self%x2=x2; self%x3=x3
    allocate(self%x2all(-1:lx2all+2),self%x3all(-1:lx3all+2))
    self%x2all=x2all; self%x3all=x3all

    self%xi_alloc_status=.true.
  end subroutine set_coords


  !> compute diffs from given grid spacing
  subroutine calc_coord_diffs(self)
    class(curvmesh) :: self
    integer :: lx1,lx2,lx3

    if (.not. self%xi_alloc_status) error stop ' attempting to compute diffs without coordinates!'

    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3    ! limits indexing verboseness, which drives me crazy

    allocate(self%dx1(0:lx1+2), self%x1i(1:lx1+1), self%dx1i(1:lx1))    
    self%dx1 = self%x1(0:lx1+2)-self%x1(-1:lx1+1)       !! computing these avoids extra message passing (could be done for other coordinates
    self%x1i(1:lx1+1) = 0.5_wp*(self%x1(0:lx1)+self%x1(1:lx1+1))
    self%dx1i=self%x1i(2:lx1+1)-self%x1i(1:lx1)

    allocate(self%dx2(0:lx2+2), self%x2i(1:lx2+1), self%dx2i(1:lx2))    
    self%dx2 = self%x2(0:lx2+2)-self%x2(-1:lx2+1)       !! computing these avoids extra message passing (could be done for other coordinates
    self%x2i(1:lx2+1) = 0.5_wp*(self%x2(0:lx2)+self%x2(1:lx2+1))
    self%dx2i=self%x2i(2:lx2+1)-self%x2i(1:lx2)

    allocate(self%dx3(0:lx3+2), self%x3i(1:lx3+1), self%dx3i(1:lx3))
    self%dx3 = self%x3(0:lx3+2)-self%x3(-1:lx3+1)
    self%x3i(1:lx3+1)=0.5_wp*(self%x3(0:lx3)+self%x3(1:lx3+1))
    self%dx3i=self%x3i(2:lx3+1)-self%x3i(1:lx3)

    self%dxi_alloc_status=.true.
  end subroutine calc_coord_diffs


  !> calculate differential lengths (units of m), needed for CFL calculations
  subroutine calc_difflengths(self)
    class(curvmesh) :: self
    real(wp), dimension(:,:,:), allocatable :: tmpdx
    integer :: lx1,lx2,lx3

    if (.not. self%dxi_alloc_status .or. .not. self%coord_alloc_status) then
      error stop ' attempting to compute differential lengths without interface diffs or metric factors!'
    end if

    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3

    allocate(tmpdx(1:lx1,1:lx2,1:lx3))
    allocate(self%dl1i(1:lx1,1:lx2,1:lx3),self%dl2i(1:lx1,1:lx2,1:lx3),self%dl3i(1:lx1,1:lx2,1:lx3))
    tmpdx=spread(spread(self%dx1i,2,lx2),3,lx3)
    self%dl1i=tmpdx*self%h1(1:lx1,1:lx2,1:lx3)
    tmpdx=spread(spread(self%dx2i,1,lx1),3,lx3)
    self%dl2i=tmpdx*self%h2(1:lx1,1:lx2,1:lx3)
    tmpdx=spread(spread(self%dx3i,1,lx1),2,lx2)
    self%dl3i=tmpdx*self%h3(1:lx1,1:lx2,1:lx3)
    deallocate(tmpdx)

    self%difflen_alloc_status=.true.
  end subroutine calc_difflengths

  !> allocate space for metric factors
  subroutine init_storage(self)
    class(curvmesh) :: self
    integer :: lx1,lx2,lx3

    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3

    if (.not. self%coord_alloc_status ) then     ! use this as a proxy for if any other coordinate-specific arrays exist
      !allocate(self%h1(1:lx1,1:lx2,1:lx3),self%h2(1:lx1,1:lx2,1:lx3),self%h3(1:lx1,1:lx2,1:lx3))
      allocate(self%h1(-1:lx1+2,-1:lx2+2,-1:lx3+2),self%h2(-1:lx1+2,-1:lx2+2,-1:lx3+2),self%h3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
      allocate(self%h1x1i(1:lx1+1,1:lx2,1:lx3),self%h2x1i(1:lx1+1,1:lx2,1:lx3),self%h3x1i(1:lx1+1,1:lx2,1:lx3))
      allocate(self%h1x2i(1:lx1,1:lx2+1,1:lx3),self%h2x2i(1:lx1,1:lx2+1,1:lx3),self%h3x2i(1:lx1,1:lx2+1,1:lx3))
      allocate(self%h1x3i(1:lx1,1:lx2,1:lx3+1),self%h2x3i(1:lx1,1:lx2,1:lx3+1),self%h3x3i(1:lx1,1:lx2,1:lx3+1))
      allocate(self%er(1:lx1,1:lx2,1:lx3,3),self%etheta(1:lx1,1:lx2,1:lx3,3),self%ephi(1:lx1,1:lx2,1:lx3,3))
      allocate(self%e1(1:lx1,1:lx2,1:lx3,3),self%e2(1:lx1,1:lx2,1:lx3,3),self%e3(1:lx1,1:lx2,1:lx3,3))
      allocate(self%Bmag(1:lx1,1:lx2,1:lx3),self%I(1:lx2,1:lx3))
      allocate(self%g1(1:lx1,1:lx2,1:lx3),self%g2(1:lx1,1:lx2,1:lx3),self%g3(1:lx1,1:lx2,1:lx3))
      allocate(self%r(1:lx1,1:lx2,1:lx3),self%theta(1:lx1,1:lx2,1:lx3),self%phi(1:lx1,1:lx2,1:lx3))
      allocate(self%alt(1:lx1,1:lx2,1:lx3),self%glon(1:lx1,1:lx2,1:lx3),self%glat(1:lx1,1:lx2,1:lx3))

      ! fixme:  there are a number of full-grid arrays that are coordinate specific to be allocated here iff we are root
      !   OR we can export that to some external procedures that collect the full grid information

      self%coord_alloc_status=.true.
    else
      error stop ' attempting to allocated space for coordinate-specific arrays when they already exist!'
    end if
  end subroutine init_storage


  !> allocate space for root-only grid quantities
  subroutine init_storage_root(self)
    class(curvmesh) :: self

    !fixme:  allocate root storage here; maybe check myid==0???  Need some protection so this does not get
    !         called from a worker...
  end subroutine init_storage_root


  subroutine calc_gridflag(self)
    class(curvmesh) :: self
    integer :: lx1,lx2,lx3

    ! error checking
    if (.not. self%geog_set_status) error stop ' attempting to compute gridflag prior to setting geographic coordinates!'

    ! sizes   
    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3

    ! grid type, assumes that each worker x1 arrays span the entire altitude range of the grid
    if (abs(self%alt(1,1,1)-self%alt(lx1,1,1))<100d3) then    !closed dipole grid
      self%gridflag=0
    else if (self%alt(1,1,1)>self%alt(2,1,1)) then            !open dipole grid with inverted structure wrt altitude
      self%gridflag=1
    else                                                      !something different (viz. non-inverted - lowest altitudes at the logical bottom of the grid)
      self%gridflag=2
    end if
  end subroutine calc_gridflag


  !> compute the number of null grid points and their indices for later use
  subroutine calc_inull(self)
    class(curvmesh) :: self
    integer :: lx1,lx2,lx3
    integer :: icount,ix1,ix2,ix3

    ! error checking, we require the geographic coords. before this is done
    if (.not. self%geog_set_status) error stop ' attempting to compute null points prior to setting geographic coordinates!'

    ! sizes   
    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3
 
    ! set null points for this simulation
    allocate(self%nullpts(1:lx1,1:lx2,1:lx3))
    self%nullpts=.false.
    where (self%alt<80e3_wp) 
      self%nullpts=.true.
    end where

    ! count needed storage for null indices
    self%lnull=0;
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          if(self%nullpts(ix1,ix2,ix3)) self%lnull=self%lnull+1
        end do
      end do
    end do
    allocate(self%inull(1:self%lnull,1:3))
    
    ! store null indices
    icount=1
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          if(self%nullpts(ix1,ix2,ix3)) then
            self%inull(icount,:)=[ix1,ix2,ix3]
            icount=icount+1
          end if
        end do
      end do
    end do

    self%null_alloc_status=.true.
  end subroutine calc_inull


  !> write the size of the grid to a file
  subroutine writesize(self,path,ID)
    class(curvmesh), intent(in) :: self
    character(*), intent(in) :: path
    integer, intent(in) :: ID
    type(hdf5_file) :: hf
    character(:), allocatable :: IDstr    ! use auto-allocation

    IDstr=' '
    write(IDstr,'(i1)') ID
    call hf%initialize(path//'simsize'//IDstr//'.h5',status='replace',action='w')
    call hf%write('/lx1',self%lx1)
    call hf%write('/lx2',self%lx2)
    call hf%write('/lx3',self%lx3)
    call hf%finalize()
  end subroutine writesize


  !> write grid coordinates (curvilinear only) to a file
  subroutine writegrid(self,path,ID)
    class(curvmesh), intent(in) :: self
    character(*), intent(in) :: path
    integer, intent(in) :: ID
    type(hdf5_file) :: hf
    character(:), allocatable :: IDstr    ! use auto-allocation

    ! at a minimum we must have allcated the coordinate arrays
    if (.not. self%xi_alloc_status) error stop ' attempting to write coordinate arrays when they have not been set yet!'

    IDstr=' '
    write(IDstr,'(i1)') ID
    call self%writesize(path,ID)
    call hf%initialize(path//'simgrid'//IDstr//'.h5',status='replace',action='w')
    call hf%write('/x1',self%x1)
    call hf%write('/x2',self%x2)
    call hf%write('/x3',self%x3)
    call hf%finalize()
  end subroutine writegrid


  subroutine writegridall(self,path,ID)
    class(curvmesh), intent(in) :: self
    character(*), intent(in) :: path
    integer, intent(in) :: ID
    type(hdf5_file) :: hf
    character(:), allocatable :: IDstr    ! use auto-allocation
    real(wp), dimension(:,:,:), allocatable :: realnullpts

    ! at a minimum we must have allcated the coordinate arrays
    if (.not. self%xi_alloc_status) error stop ' attempting to write coordinate arrays when they have not been set yet!'

    IDstr=' '
    write(IDstr,'(i1)') ID
    call self%writesize(path,ID)
    call hf%initialize(path//'simgrid'//IDstr//'.h5',status='replace',action='w')
    call hf%write('/x1',self%x1)
    call hf%write('/x1i',self%x1i)
    call hf%write('/dx1b',self%dx1)
    call hf%write('/dx1h',self%dx1i)

    call hf%write('/x2',self%x2)
    call hf%write('/x2i',self%x2i)
    call hf%write('/dx2b',self%dx2)
    call hf%write('/dx2h',self%dx2i)

    call hf%write('/x3',self%x3)
    call hf%write('/x3i',self%x3i)
    call hf%write('/dx3b',self%dx3)
    call hf%write('/dx3h',self%dx3i)

    call hf%write('/h1',self%h1)
    call hf%write('/h2',self%h2)
    call hf%write('/h3',self%h3)

    call hf%write('/h1x1i',self%h1x1i)
    call hf%write('/h2x1i',self%h2x1i)
    call hf%write('/h3x1i',self%h3x1i)

    call hf%write('/h1x2i',self%h1x2i)
    call hf%write('/h2x2i',self%h2x2i)
    call hf%write('/h3x2i',self%h3x2i)

    call hf%write('/h1x3i',self%h1x3i)
    call hf%write('/h2x3i',self%h2x3i)
    call hf%write('/h3x3i',self%h3x3i)

    call hf%write('/gx1',self%g1)
    call hf%write('/gx2',self%g2)
    call hf%write('/gx3',self%g3)

    call hf%write('/alt',self%alt)
    call hf%write('/glon',self%glon)
    call hf%write('/glat',self%glat)

    call hf%write('/Bmag',self%Bmag)
    call hf%write('/I',self%I)
    allocate(realnullpts(1:self%lx1,1:self%lx2,1:self%lx3))
    realnullpts=0.0
    where (self%nullpts)
      realnullpts=1._wp
    end where
    call hf%write('/nullpts',realnullpts)
    deallocate(realnullpts)

    call hf%write('/e1',self%e1)
    call hf%write('/e2',self%e2)
    call hf%write('/e3',self%e3)

    call hf%write('/er',self%er)
    call hf%write('/etheta',self%etheta)
    call hf%write('/ephi',self%ephi)

    call hf%write('/r',self%r)
    call hf%write('/theta',self%theta)
    call hf%write('/phi',self%phi)

    call hf%finalize()

  end subroutine writegridall


  !> type destructor; written generally, viz. as if it is possible some grid pieces are allocated an others are not
  subroutine destructor(self)
    type(curvmesh) :: self

    ! deallocation statements here; always check allocation status flags first...
    if (self%xi_alloc_status) deallocate(self%x1,self%x2,self%x3,self%x2all,self%x3all)    ! these are from set_coords
    if (self%dxi_alloc_status) then                                  ! from calc_coord_diffs
      deallocate(self%dx1,self%x1i,self%dx1i)
      deallocate(self%dx2,self%x2i,self%dx2i)
      deallocate(self%dx3,self%x3i,self%dx3i)
    end if
    if (self%difflen_alloc_status) deallocate(self%dl1i,self%dl2i,self%dl3i)    ! from calc_difflengths

    ! coordinate-specific arrays set by type extensions
    if (self%coord_alloc_status) then
      deallocate(self%h1,self%h2,self%h3,self%er,self%etheta,self%ephi,self%e1,self%e2,self%e3)
      deallocate(self%r,self%theta,self%phi)
      deallocate(self%h1x1i,self%h2x1i,self%h3x1i)
      deallocate(self%h1x2i,self%h2x2i,self%h3x2i)
      deallocate(self%g1,self%g2,self%g3)
      deallocate(self%Bmag,self%I)
      deallocate(self%alt,self%glon,self%glat)
    end if

    if (self%null_alloc_status) then
      deallocate(self%nullpts,self%inull)
    end if

    ! let the user know that the destructor indeed ran
    print*, '  curvmesh destructor completed successfully'
  end subroutine destructor

end module meshobj
