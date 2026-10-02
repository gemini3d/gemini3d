module audit_driver_fixture
use inputdataobj, only: inputdata
use phys_consts, only: wp
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use timeutils, only: shift_datetime,elapsed_seconds
implicit none
type, extends(inputdata) :: linear_driver
  integer :: origin(3),loads=0
  real(wp) :: origin_ut,value
contains
  procedure :: init=>fixture_init
  procedure :: set_coordsi=>fixture_coords
  procedure :: load_size=>fixture_noop
  procedure :: load_grid=>fixture_noop
  procedure :: load_data=>fixture_load
  procedure :: spaceinterp=>fixture_space
end type
contains
subroutine fixture_init(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
  class(linear_driver),intent(inout)::self
  type(gemini_cfg),intent(in)::cfg
  character(*),intent(in)::sourcedir
  class(curvmesh),intent(in)::x
  real(wp),intent(in)::dtmodel,dtdata,UTsec
  integer,intent(in)::ymd(3)
  self%origin=cfg%ymd0;self%origin_ut=cfg%UTsec0
  self%l0D=1;self%l1Dax1=0;self%l1Dax2=0;self%l1Dax3=0
  self%l2Dax12=0;self%l2Dax13=0;self%l2Dax23=0;self%l3D=0
  self%lc1i=1;self%lc2i=1;self%lc3i=1
  allocate(self%data0Di(1,2),self%data0Dinow(1))
  self%data0Di=0;self%flagalloc=.true.;self%flagdoinput=.true.
  call self%set_cadence(dtdata)
  call self%set_source(sourcedir)
  call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
end subroutine
subroutine fixture_coords(self,cfg,x)
  class(linear_driver),intent(inout)::self
  type(gemini_cfg),intent(in)::cfg
  class(curvmesh),intent(in)::x
  self%flagcoordsi=.true.
end subroutine
subroutine fixture_noop(self)
  class(linear_driver),intent(inout)::self
end subroutine
subroutine fixture_load(self,t,dtmodel,ymdtmp,UTsectmp)
  class(linear_driver),intent(inout)::self
  real(wp),intent(in)::t,dtmodel
  integer,intent(inout)::ymdtmp(3)
  real(wp),intent(inout)::UTsectmp
  ymdtmp=self%ymdref(:,2);UTsectmp=self%UTsecref(2)
  call shift_datetime(self%dt,ymdtmp,UTsectmp)
  self%value=5+2*elapsed_seconds(self%origin,self%origin_ut,ymdtmp,UTsectmp)
  self%loads=self%loads+1
end subroutine
subroutine fixture_space(self)
  class(linear_driver),intent(inout)::self
  self%data0Di(:,1)=self%data0Di(:,2)
  self%data0Di(:,2)=self%value
end subroutine
end module
