program qualification_coverage
use audit_driver_fixture
use interpolation, only: coverage_mask,interp1
use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
implicit none
type(linear_driver) :: data
character(20) :: mode
real(wp) :: result(4),values(2)
call get_command_argument(1,mode)
data%l1Dax1=1; data%l1Dax2=0; data%l1Dax3=0
data%l2Dax12=0; data%l2Dax13=0; data%l2Dax23=0; data%l3D=0
allocate(data%coord1(2),data%coord1iax1(4))
data%coord1=[0._wp,1._wp]
data%coord1iax1=[0._wp,0.5_wp,1._wp,1._wp]
values=0 ! physical zero must be distinguishable from absent coverage
if (trim(mode)=='nan') then
  values(1)=ieee_value(0._wp,ieee_quiet_nan)
  result=interp1(data%coord1,values,data%coord1iax1)
  stop ! acceptance here must make the WILL_FAIL test fail
endif
if (trim(mode)=='flagged'.or.trim(mode)=='reject') data%coord1iax1(4)=2._wp
data%allow_missing_spatial=trim(mode)=='flagged'
call data%validate_spatial_coverage()
if (trim(mode)=='reject') stop
if (trim(mode)=='flagged') then
  if (data%spatial_missing_count/=1.or.data%coverage(1)%valid(4)) error stop 'missing flag lost'
  if (.not.all(data%coverage(1)%valid(1:3))) error stop 'physical zero marked missing'
else
  if (data%spatial_missing_count/=0) error stop 'covered sites marked missing'
endif
if (.not.all(coverage_mask([5._wp],[-100._wp,100._wp]))) error stop 'singleton invariant axis'
deallocate(data%coord1,data%coord1iax1)
print *, 'Coverage policy passed: ',trim(mode)
end program
