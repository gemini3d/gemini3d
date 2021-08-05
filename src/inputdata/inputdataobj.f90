module inputdataobj

use phys_consts, only : wp
use meshobj, only : curvmesh
use interpolation, only : interp1,interp2,interp3

implicit none (type, external)
public

!> this is a generic class for an data object being input into the model and interpolated in space and time
type, abstract :: inputdata
  !! here we store data that have already been received but not yet interpolated
  real(wp), dimension(:), pointer :: coord1,coord2,coord3     ! coordinates for the source data (interpolant coords)
  integer :: lc1,lc2,lc3                                      ! dataset length along the 3 coordinate axes
  real(wp), dimension(:), pointer :: data0D
  real(wp), dimension(:,:), pointer :: data1Dax1,data1Dax2,data1Dax3
  real(wp), dimension(:,:), pointer :: data2Dax23,data2Dax12,data2Dax13
  real(wp), dimension(:,:), pointer :: data3D


  !! here we store data that have already been spatially interpolated
  real(wp), dimension(:,:), pointer :: data0Di                    ! array for storing a "stack" of scalar data (only interpolated in time)
                                                                     !  last axis is for prev,next copies of data for interp in time
                                                                     !  second to last axis is for number of datasets of this dimension
  integer :: l0D                                                     ! length/number of scalar datasets
  real(wp), dimension(:,:,:), pointer :: data1Dax1i                  ! array for storing series of 1D data, array varies along non-singleton axis
  real(wp), dimension(:,:,:), pointer :: data1Dax2i,data1Dax3i       ! 1D data arrays varying along coordinates (axes) 2 and 3 
  integer :: l1Dax1,l1Dax2,l1Dax3
  real(wp), dimension(:,:,:,:), pointer :: data2Dax23i               ! array for storing series of 2D data, varies along two non-singleton axes
  real(wp), dimension(:,:,:,:), pointer :: data2Dax12i,data2dax13i   !2D arrays varying along 1,2 and 1,3 axes 
  integer :: l2Dax23,l2Dax12,l2Dax13
  real(wp), dimension(:,:,:,:,:), pointer :: data3Di                 ! array for storing series of 3D data
  integer :: lparm3D
  real(wp), dimension(:,:,:), pointer :: coord1i,coord2i,coord3i     ! coordinates of the interpolation sites
  integer :: lc1i,lc2i,lc3i                                          ! dataset length along the 3 coordinate axes

  real(wp), dimension(:), pointer :: data0Dinow 
  real(wp), dimension(:,:), pointer :: data1Dax1inow
  real(wp), dimension(:,:), pointer :: data1Dax2inow,data1Dax3inow
  real(wp), dimension(:,:,:), pointer :: data2Dax23inow
  real(wp), dimension(:,:,:), pointer :: data2Dax12inow,data2dax13inow  
  real(wp), dimension(:,:,:,:), pointer :: data3Dinow

  real(wp), dimension(2) :: tref                                     ! times for two input frames bracketting current time
  real(wp) :: tnow                                                   ! time corresponding to data in *now arrays, viz current time insofar as this object knows

  contains
    procedure :: setsizes              ! initiate sizes for coordinate axes and number of datasets of different dimensionality
    procedure :: setcoords             ! fill interpolant coordinate arrays
    procedure :: init_storage          ! wrapper routine to set up arrays once sizes are known/set
    procedure :: update                ! check to see if new file needs to be read and read accordingly (will need to call deferred loaddata)
    procedure :: timeinterp            ! interpolate in time based on data presently loaded into spatial arrays
    procedure :: spaceinterp           ! interpolate spatially
    procedure :: dissociate_pointers   ! clear out memory and reset and allocation status flags

    procedure(initproc), deferred :: init                 ! create storage (sub), and prime data as needed
    procedure(coordsetproc), deferred :: setcoordsi       ! use grid data to compute coordinates of the interpolation sites
    procedure(loadproc), deferred :: loaddata             ! read data from file (possibly one array at a time) and spatially interpolate and store it in the appropriate arrays
end type inputdata

!> interfaces for deferred procedures
abstract interface
  subroutine initproc(self)
    import inputdata
    class(inputdata), intent(inout) :: self
  end subroutine init
  subroutine (coordsetproc(self,x)
    import inputdata
    class(inputdata), intent(inout) :: self
    class(curvmesh), intent(in)  :: x
  end subroutine coordsetproc
  subroutine loadproc(self)
    import inputdata
    class(inputdata), intent(inout) :: self
  end subroutine loadproc
end interface

contains
  !> use data stored in input arrays to interpolate onto grid sites for "next" dataset
  subroutine spaceinterp(self)
    class(inputdata),intent(inout) :: self
    integer :: iparm
    real(wp), dimension(:), allocatable :: tempdata

    !> 1D arrays varying along the 1-axis
    if (associated(data1Dax1i)) then
      allocate(tempdata(self%lc1i))
      do iparm=1,l1Dax1
        tempdata(:)=interp1(coord1,data1Dax1(:,iparm),coord1i)
        self%data1Dax1i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 1D arrays varying along the 2-axis
    if (associated(data1Dax2i)) then
      allocate(tempdata(self%lc2i))
      do iparm=1,l1Dax2
        tempdata(:)=interp1(coord2,data1Dax2(:,iparm),coord2i)
        self%data1Dax2i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 1D arrays varying along the 3-axis
    if (associated(data1Dax3i)) then
      allocate(tempdata(self%lc3i))
      do iparm=1,l1Dax3
        tempdata(:)=interp1(coord3,data1Dax3(:,iparm),coord3i)
        self%data1Dax3i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 2D arrays varying along the 2,3 axes
    if (associated(data2Dax23i)) then
      allocate(tempdata(self%lc2i*self%lc3i))
      do iparm=1,l2D23
        tempdata(:)=interp2(coord2,coord3,data2Dax23(:,:,iparm),coord2i,coord3i)
        self%data2Dax23i(:,:,iparm,2)=reshape(tempdata,[lc2i,lc3i])
      end do
      deallocate(tempdata)
    end if

    !> 2D arrays varying along the 1,2 axes
    if (associated(data2Dax12i)) then
      allocate(tempdata(self%lc1i*self%lc2i))
      do iparm=1,l2D12
        tempdata(:)=interp2(coord1,coord2,data2Dax12(:,:,iparm),coord1i,coord2i)
        self%data2Dax12i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i])
      end do
      deallocate(tempdata)
    end if

    !> 2D arrays varying along the 1,3 axes
    if (associated(data2Dax13i)) then
      allocate(tempdata(self%lc1i*self%lc3i))
      do iparm=1,l2D13
        tempdata(:)=interp2(coord1,coord3,data2Dax13(:,:,iparm),coord1i,coord3i)
        self%data2Dax13i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc3i])
      end do
      deallocate(tempdata)
    end if

    !> 3D arrays varying along all axes (obv.)
    if (associated(data3Di)) then
      allocate(tempdata(self%lc1i*self%lc2i*self%lc3i))
      do iparm=1,l3D
        tempdata(:)=interp3(coord1,coord2,coord3,data3D(:,:,:,iparm),coord1i,coord2i,coord3i)
        self%data3Di(:,:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i,lc3i])
      end do
      deallocate(tempdata)
    end if 
  end subroutine spaceinterp


  !> interpolate data in time
  subroutine timeinterp(self,t,dt)
    class(inputdata), intent(inout) :: self
    real(wp), intent(in) :: t
    real(wp) :: slope
    integer :: ic1,ic2,ic3

    !> interpolate scalars in time
    if (associated(data0Di)) then
      slope=(self%data0Di(iparm,2)-self%data0Di(iparm,1))/(self%tref(2)-self%tref(1))
      self%data0Dinow(ic1,iparm)=self%data0Di(iparm,1)+slope*(t+dt/2 -self%tref(1))
    end if

    !> 1D arrays varying along the 1-axis
    if (associated(data1Dax1i)) then
      do iparm=1,l1Dax1
        do ic1=1,lc1i
          slope=(self%data1Dax1i(ic1,iparm,2)-self%data1Dax1i(ic1,iparm,1))/(self%tref(2)-self%tref(1))
          self%data1Dax1inow(ic1,iparm)=self%data1Dax1i(ic1,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end if

    !> 1D arrays varying along the 2-axis
    if (associated(data1Dax2i)) then
      do iparm=1,l1Dax2
        do ic2=1,lc2i
          slope=(self%data1Dax2i(ic2,iparm,2)-self%data1Dax2i(ic2,iparm,1))/(self%tref(2)-self%tref(1))
          self%data1Dax2inow(ic2,iparm)=self%data1Dax2i(ic2,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end if

    !> 1D arrays varying along the 3-axis
    if (associated(data1Dax3i)) then
      do iparm=1,l1Dax3
        do ic3=1,lc3i
          slope=(self%data1Dax3i(ic3,iparm,2)-self%data1Dax3i(ic3,iparm,1))/(self%tref(2)-self%tref(1))
          self%data1Dax3inow(ic3,iparm)=self%data1Dax3i(ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end if

    !> 2D arrays varying along the 2,3 axes
    if (associated(data2Dax23i)) then
      do iparm=1,l2D23
        do ic3=1,lc3i
          do ic2=1,lc2i
            slope=(self%data2Dax23i(ic2,ic3,iparm,2)-self%data2Dax23i(ic2,ic3,iparm,1))/(self%tref(2)-self%tref(1))
            self%data2Dax23inow(ic2,ic3,iparm)=self%data2Dax23i(ic2,ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
          end do
        end do
      end do
    end if

    !> 2D arrays varying along the 1,2 axes
    if (associated(data2Dax12i)) then
      do iparm=1,l2D12
        do ic2=1,lc2i
          do ic1=1,lc1i
            slope=(self%data2Dax12i(ic1,ic2,iparm,2)-self%data2Dax12i(ic1,ic2,iparm,1))/(self%tref(2)-self%tref(1))
            self%data2Dax12inow(ic1,ic2,iparm)=self%data2Dax12i(ic1,ic2,iparm,1)+slope*(t+dt/2 -self%tref(1))
          end do
        end do
      end do
    end if

    !> 2D arrays varying along the 1,3 axes
    if (associated(data2Dax13i)) then
      do iparm=1,l2D13
        do ic3=1,lc3i
          do ic1=1,lc1i
            slope=(self%data2Dax13i(ic1,ic3,iparm,2)-self%data2Dax13i(ic1,ic3,iparm,1))/(self%tref(2)-self%tref(1))
            self%data2Dax13inow(ic1,ic3,iparm)=self%data2Dax13i(ic1,ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
          end do
        end do
      end do
    end if

    !> 3D arrays varying along all axes (obv.)
    if (associated(data3Di)) then
      do iparm=1,l3D
        do ic3=1,lc3i
          do ic2=1,lc2i
            do ic1=1,lc1i
              slope=(self%data3Di(ic1,ic2,ic3,iparm,2)-self%data3Di(ic1,ic2,ic3,iparm,1))/(self%tref(2)-self%tref(1))
              self%data3Dinow(ic1,ic2,ic3,iparm)=self%data3Di(ic1,ic2,ic3,iparm,1)+slope*(t+dt/2-self%tref(1))
            end do
          end do
        end do
      end do
    end if 

    !> update current time in object
    self%tnow=t+dt/2
  end subroutine timeinterp
end module inputdataobj
