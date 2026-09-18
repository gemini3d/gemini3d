
program audit_driver
use audit_driver_fixture
use meshobj_cart, only: cartmesh
implicit none
type(linear_driver)::driver
type(gemini_cfg)::cfg
type(cartmesh)::x
cfg%ymd0=[2019,12,31];cfg%UTsec0=86370
! Restart two minutes after the origin, across the year boundary.
call driver%init(cfg,'synthetic',x,2._wp,60._wp,[2020,1,1],90._wp)
if(any(abs(driver%tref-[120._wp,180._wp])>1e-12_wp))error stop 'midnight priming interval'
if(driver%loads/=2)error stop 'priming must load exactly two frames'
call driver%update(cfg,2._wp,120._wp,x,[2020,1,1],90._wp)
if(abs(driver%data0Dinow(1)-247._wp)>1e-12_wp)error stop 'linear driver interpolation'
! Skip multiple driver cadences; catch up before interpolating.
call driver%update(cfg,20._wp,370._wp,x,[2020,1,1],340._wp)
if(abs(driver%data0Dinow(1)-765._wp)>1e-12_wp)error stop 'multiple cadence catch-up'
if(any(driver%ymdref(:,2)/=[2020,1,1]))error stop 'driver date mismatch'
if(abs(driver%UTsecref(2)-390._wp)>1e-12_wp)error stop 'driver clock mismatch'
deallocate(driver%data0Di,driver%data0Dinow)
driver%flagalloc=.false.
print *, 'Cross-midnight driver priming and catch-up passed'
end program
