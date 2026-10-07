module neutral_abi_fixture
use, intrinsic :: iso_c_binding, only: c_ptr,c_loc,c_f_pointer,c_int
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp
use gemini3d, only: gemini_work
use gemini3d_config, only: gemini_cfg
use meshobj_cart, only: cartmesh
use grid, only: set_total_grid_sizes,set_subgrid_sizes
use neutral, only: neutral_info_alloc,neutral_info_dealloc
implicit none
type(gemini_cfg), target, save :: cfg
type(cartmesh), target, save :: x
type(gemini_work), target, save :: work
contains
subroutine fixture(cfgC,xC,workC,version) bind(C,name="audit_neutral_fixture")
type(c_ptr), intent(out) :: cfgC,xC,workC
integer(c_int), intent(in) :: version
call set_total_grid_sizes(2,1,1)
call set_subgrid_sizes(1,1)
cfg%activ=[100._wp,100._wp,4._wp]
cfg%msis_version=version
x%lx1=2; x%lx2=1; x%lx3=1; x%lnull=0
allocate(x%glat(2,1,1),x%glon(2,1,1),x%alt(2,1,1),work%atmos)
x%glat=65; x%glon=210; x%alt(:,1,1)=[200e3_wp,300e3_wp]
call neutral_info_alloc(work%atmos)
work%atmos%flagprojections=.true.
work%atmos%proj_ealt_e1=1; work%atmos%proj_eglat_e1=0; work%atmos%proj_eglon_e1=0
work%atmos%proj_ealt_e2=0; work%atmos%proj_eglat_e2=1; work%atmos%proj_eglon_e2=0
work%atmos%proj_ealt_e3=0; work%atmos%proj_eglat_e3=0; work%atmos%proj_eglon_e3=1
cfgC=c_loc(cfg); xC=c_loc(x); workC=c_loc(work)
end subroutine

subroutine check(workC,v2,v3) bind(C,name="audit_neutral_check")
type(c_ptr), intent(in) :: workC
real(wp), intent(in) :: v2,v3
type(gemini_work), pointer :: w
call c_f_pointer(workC,w)
if (.not.all(ieee_is_finite(w%atmos%nn)) .or. .not.all(ieee_is_finite(w%atmos%Tn))) &
  error stop "nonfinite empirical atmosphere through C ABI"
if (any(w%atmos%nn<=0) .or. any(w%atmos%Tn<100)) error stop "C ABI did not compute empirical atmosphere"
if (any(w%atmos%nn/=w%atmos%nnBG) .or. any(w%atmos%Tn/=w%atmos%TnBG)) error stop "background not aggregated"
if (any(abs(w%atmos%vn1-w%atmos%vn1BG)>1e-10_wp)) error stop "wind 1 not aggregated"
if (any(abs(w%atmos%vn2-(w%atmos%vn2BG-v2))>1e-10_wp)) error stop "wind 2 drift not applied exactly once"
if (any(abs(w%atmos%vn3-(w%atmos%vn3BG-v3))>1e-10_wp)) error stop "wind 3 drift not applied exactly once"
end subroutine

subroutine release(workC) bind(C,name="audit_neutral_release")
type(c_ptr), intent(in) :: workC
type(gemini_work), pointer :: w
call c_f_pointer(workC,w)
call neutral_info_dealloc(w%atmos)
deallocate(w%atmos,x%glat,x%glon,x%alt)
end subroutine
end module
