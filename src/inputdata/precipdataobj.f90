module precipdataobj

use phys_consts, only: wp
use inputdataobj, only: inputdata
use meshobj, only: curvmesh
use config, only: gemini_cfg
use io, only: ...
implicit none (type, external)

type, extends(inputdata) :: precipdata
  ! coordinate for input precipitation data, and storage
  real(wp), dimension(:), pointer :: mlonp,mlatp
  integer, pointer :: llon,llat
  real(wp), dimension(:,:), pointer :: Qp,E0p

  contains
    final :: destructor
    procedure(initproc) :: init=>init_precip
    procedure(coordsisetproc) :: set_coordsi=>set_coordsi_precip
    procedure(loadproc) :: load_data=>load_data_precip
    procedure(loadgridproc) :: load_grid=>load_grid_precip
end type precipdata

contains
  !> set pointers to appropriate data arrays (taking into account dimensionality of the problem)
  subroutine init_precip(self,cfg,sourcedir,x,dtdata,ymd,UTsec)
    class(precipdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg                 ! gemini config type for optional params
    character, dimension(:), intent(in) :: sourcedir    ! directory for precipitation data input
    class(curvmesh), intent(in) :: x                    ! curvmesh object
    real(wp), intent(in) :: dtdata                      ! cadence for input data from config.nml
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation

    ! read the simulation size from the source directory and allocate arrays
    call get_simsize()
    call self%set_sizes(lc1,lc2,lc3, &
                       0, &
                       0,0,0 &
                       2,0,0 &
                       0, &
                       x )
    call self%init_storage()
    call self%set_cadence(dtdata)
    call self%set_name('electron precipitation')

    ! set local pointers grid pointers and assign input data grid
    self%mlonp=>self%coord2; self%mlatp=>self%coord3;
    llon=>lc2; llat=>lc3
    call self%load_grid_precip(mlonp,mlatp)

    ! set input data array pointers to faciliate easy to read input code
    Qp=>data2Dax23(:,:,1)
    E0p=>data2Dax23(:,:,2)
    Qpiprev=>data2Dax23i(:,:,1,1)
    Qpinext=>data2Dax23i(:,:,1,2)
    E0iprev=>data2Dax23i(:,:,2,1)
    E0inext=>data2Dax23i(:,:,2,2)
    Qpinow=>data2Dax23inow(:,:,1)
    E0inow=>data2Dax23inow(:,:,2)

    ! prime input data
    call self%prime_data(cfg,x,ymd,UTsec)
  end subroutine init_precip
end module precipdataobj
