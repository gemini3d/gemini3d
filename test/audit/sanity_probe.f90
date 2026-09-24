! Audit addition 2026-09-16. Apache-2.0.
program audit_sanity
use phys_consts, only: wp
use sanity_check, only: check_finite_plasma, check_finite_output
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
implicit none
real(wp) :: ns(6,6,6,7), v1(6,6,6,7), v2(6,6,6,7), v3(6,6,6,7), ts(6,6,6,7)
real(wp) :: phi(6,6,6), j(6,6,6), nan
character(128) :: field, species_text, kind
character(2048) :: outdir
integer :: species
call get_command_argument(1,field); call get_command_argument(2,species_text)
call get_command_argument(3,kind); call get_command_argument(4,outdir)
read(species_text,*) species
nan=ieee_value(0._wp,ieee_quiet_nan)
ns=nan; v1=nan; ts=nan; phi=nan
ns(3:4,3:4,3:4,:)=1e11_wp; v1(3:4,3:4,3:4,:)=0; ts(3:4,3:4,3:4,:)=1000
phi(3:4,3:4,3:4)=0; j=phi; v2=v1; v3=v1
if (trim(kind)=='null_partition') ns(3:4,3:4,3:4,:)=1e-20_wp
if (trim(kind)=='nan') then
 select case(trim(field))
 case('ns'); ns(3,3,3,species)=nan
 case('v1'); v1(3,3,3,species)=nan
 case('v2'); v2(3,3,3,species)=nan
 case('v3'); v3(3,3,3,species)=nan
 case('ts'); ts(3,3,3,species)=nan
 end select
elseif (trim(kind)=='negative') then
 if(trim(field)=='ts')then
  ts(3,3,3,species)=-1
 else
  ns(3,3,3,species)=-1
 endif
endif
if (trim(field)=='v2' .or. trim(field)=='v3') then
 call check_finite_output(trim(outdir),0._wp,0,v2,v3,ns,v1,ts,phi,j,j,j)
else
 call check_finite_plasma(trim(outdir),ns,v1,ts)
endif
print *, "ACCEPTED"
end program
