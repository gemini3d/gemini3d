program qualification_capacitance_mode
use phys_consts, only: wp
use gemini3d_config, only: gemini_cfg
use meshobj_cart, only: cartmesh
use potential_comm, only: validate_capacitance_mode
use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
implicit none
type(gemini_cfg) :: cfg
type(cartmesh) :: x
character(20) :: mode
integer :: boundary
call get_command_argument(1,mode)
cfg%flagcap=1;cfg%potsolve=1;boundary=0
x%lx1=3;x%lx2=3;x%lx3=3;x%lx2all=3;x%lx3all=3
allocate(x%h1(3,3,3),x%h2(3,3,3),x%h3(3,3,3))
x%h1=1;x%h2=1;x%h3=1
select case(trim(mode))
case('supported')
case('2d');x%lx2all=1
case('field_resolved');cfg%potsolve=3
case('dirichlet');boundary=1
case('current_bc');boundary=2
case('curved');x%h2=1.2_wp
case('nan');x%h3(1,1,1)=ieee_value(1._wp,ieee_quiet_nan)
case default;error stop 'Unknown test mode'
end select
call validate_capacitance_mode(cfg,x,boundary)
deallocate(x%h1,x%h2,x%h3)
! Negative modes return success if a regression incorrectly permits them;
! CTest WILL_FAIL therefore reports failure rather than a false pass.
end program
