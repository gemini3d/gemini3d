program neutral_background_test
use phys_consts, only: wp
use grid, only: lx1,lx2,lx3,set_total_grid_sizes,set_subgrid_sizes
use meshobj_cart, only: cartmesh
use gemini3d_config, only: gemini_cfg
use neutral, only: neutral_info,neutral_info_alloc,neutral_info_dealloc
use neutral_background, only: neutral_background_fileinput
use neutraldataBGobj, only: neutraldataBG
implicit none
type(cartmesh) :: x
type(gemini_cfg) :: cfg
type(neutral_info) :: atmos
type(neutraldataBG), allocatable :: bg
character(40) :: mode
real(wp), allocatable :: expected_n(:),expected_t(:),expected_v(:)
integer :: i

call get_command_argument(1,mode)
call set_total_grid_sizes(5,1,1)
if (mode=='singleton' .or. mode=='singleton_uncovered') call set_total_grid_sizes(1,1,1)
call set_subgrid_sizes(1,1)
x%lx1=lx1; x%lx2=lx2; x%lx3=lx3
allocate(x%alt(lx1,1,1),bg)
allocate(bg%data3Dinow(lx1,1,1,9))
bg%natminow=>bg%data3Dinow
bg%flagdoinput=.false. ! exercise copyout using a prepared, already interpolated frame
bg%altpmax=300e3_wp
call neutral_info_alloc(atmos)
atmos%flagprojections=.true.
atmos%proj_ealt_e1=1; atmos%proj_eglat_e1=0; atmos%proj_eglon_e1=0
atmos%proj_ealt_e2=0; atmos%proj_eglat_e2=1; atmos%proj_eglon_e2=0
atmos%proj_ealt_e3=0; atmos%proj_eglat_e3=0; atmos%proj_eglon_e3=1
atmos%vn1=-999; atmos%vn2=-999; atmos%vn3=-999
allocate(expected_n(lx1),expected_t(lx1),expected_v(lx1))
do i=1,lx1
  x%alt(i,1,1)=100e3_wp*i
  bg%natminow(i,1,1,1:5)=16._wp/2**(i-1)
  bg%natminow(i,1,1,6:8)=10._wp*i
  bg%natminow(i,1,1,9)=500._wp+10*i
end do
expected_n=bg%natminow(:,1,1,1)
expected_t=bg%natminow(:,1,1,9)
expected_v=bg%natminow(:,1,1,6)

select case (mode)
case ('covered','endpoint','singleton')
  bg%altpmax=maxval(x%alt)
  if (mode=='covered') bg%altpmax=bg%altpmax+1
case ('descending_covered','descending')
  x%alt(:,1,1)=x%alt(lx1:1:-1,1,1)
  bg%natminow(:,1,1,:)=bg%natminow(lx1:1:-1,1,1,:)
  expected_n=expected_n(lx1:1:-1)
  expected_t=expected_t(lx1:1:-1)
  expected_v=expected_v(lx1:1:-1)
  if (mode=='descending_covered') then
    bg%altpmax=600e3_wp
  else
    expected_t(1:2)=expected_t(3)
    expected_v(1:2)=expected_v(3)
  end if
case ('closed','closed_covered')
  x%alt(:,1,1)=[100e3_wp,200e3_wp,300e3_wp,200e3_wp,100e3_wp]
  bg%natminow(4:5,1,1,1:5)=bg%natminow(2:1:-1,1,1,1:5)
  expected_n=[16._wp,8._wp,4._wp,8._wp,16._wp]
  if (mode=='closed') then
    bg%altpmax=200e3_wp
    ! Preserve the existing maximum-x1 hemisphere convention at the apex.
    expected_v(3)=expected_v(4); expected_t(3)=expected_t(4)
  end if
case ('ascending','zero','belowground')
  expected_t(4:5)=expected_t(3)
  expected_v(4:5)=expected_v(3)
  if (mode=='zero') then
    bg%natminow(:,1,1,4)=0
  else if (mode=='belowground') then
    x%alt(1,1,1)=-10e3_wp
    expected_n(1)=expected_n(2); expected_t(1)=expected_t(2)
  end if
case ('lower_coverage')
  bg%altpmax=100e3_wp
case ('upper_coverage')
  x%alt(:,1,1)=x%alt(lx1:1:-1,1,1)
  bg%altpmax=100e3_wp
case ('zero_denominator')
  bg%natminow(2,1,1,4)=0
case ('negative_density')
  bg%natminow(2,1,1,4)=-1
case ('underground')
  x%alt=-1
case ('temperature')
  bg%natminow(3,1,1,9)=0
case ('singleton_uncovered')
  bg%altpmax=50e3_wp
case default
  error stop 'unknown neutral test mode'
end select

call neutral_background_fileinput(1._wp,0._wp,cfg,[2020,1,1],0._wp,x,atmos,bg)
if (any(abs(atmos%nnBG(:,1,1,1)-expected_n)>1e-12_wp)) error stop 'density extrapolation'
if (any(abs(atmos%TnBG(:,1,1)-expected_t)>1e-12_wp)) error stop 'temperature extrapolation'
if (any(abs(atmos%vn2BG(:,1,1)-expected_v)>1e-12_wp)) error stop 'fresh wind 2 extrapolation'
if (any(abs(atmos%vn3BG(:,1,1)-expected_v)>1e-12_wp)) error stop 'fresh wind 3 extrapolation'
expected_v=expected_v*(0.5_wp+0.5_wp*tanh((x%alt(:,1,1)-150e3_wp)/10e3_wp))
if (any(abs(atmos%vn1BG(:,1,1)-expected_v)>1e-12_wp)) error stop 'fresh wind 1 extrapolation'
if (mode=='zero') then
  if (any(atmos%nnBG(:,1,1,4)/=0)) error stop 'absent species must remain absent'
end if
if (any(atmos%vn1/=-999) .or. any(atmos%vn2/=-999) .or. any(atmos%vn3/=-999)) &
  error stop 'background update changed aggregate winds'
call neutral_info_dealloc(atmos)
deallocate(x%alt,bg)
print *, 'neutral background passed: ',trim(mode)
end program neutral_background_test
