program audit_cfl
use phys_consts, only: wp,lsp
use temporal, only: cflcalc
use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
implicit none
real(wp) :: ts(-1:4,-1:3,-1:3,lsp),v1(-1:4,-1:3,-1:3,lsp),v2(-1:4,-1:3,-1:3,lsp),v3(-1:4,-1:3,-1:3,lsp)
real(wp) :: dl1(2,1,1),dl2(2,1,1),dl3(2,1,1),cfl,dt
character(20) :: mode
call get_command_argument(1,mode)
ts=300;v1=0;v2=0;v3=0;dl1=1e6_wp;dl2=1;dl3=1;dt=0.1_wp
v2(2,1,1,7)=7
if(trim(mode)=='nan')v2(2,1,1,7)=ieee_value(0._wp,ieee_quiet_nan)
if(trim(mode)=='metric')dl2=0
if(trim(mode)=='temperature')ts(1,1,1,1)=-1
if(trim(mode)=='step')dt=0
call cflcalc(ts,v1,v2,v3,dl1,dl2,dl3,dt,cfl)
if(abs(cfl-0.7_wp)>1e-12_wp)error stop 'CFL advection oracle mismatch'
v2(2,1,1,7)=20
call cflcalc(ts,v1,v2,v3,dl1,dl2,dl3,dt,cfl)
if(cfl<=1)error stop 'updated drift did not invalidate old step'
print *,'CFL checks passed'
end program
